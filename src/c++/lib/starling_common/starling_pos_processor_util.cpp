//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

/// \file
/// \author Chris Saunders
///


#include "CandidateAlignment.hh"
#include "normalizeAlignment.hh"
#include "starling_pos_processor_indel_util.hh"
#include "starling_read_filter_shared.hh"
#include "starling_read_util.hh"

#include "blt_util/align_path.hh"
#include "blt_util/log.hh"
#include "blt_util/qscore.hh"
#include "common/Exceptions.hh"
#include "htsapi/bam_header_util.hh"
#include "htsapi/align_path_bam_util.hh"
#include "htsapi/bam_seq.hh"
#include "starling_common/starling_pos_processor_util.hh"

#include <cassert>

#include <sstream>



std::vector<std::reference_wrapper<const bam_hdr_t> >
registerAlignments(
    const std::vector<std::string>& alignmentFilename,
    const std::vector<unsigned>& registrationIndices,
    HtsMergeStreamer& streamData)
{
    const unsigned alignmentFileCount(alignmentFilename.size());
    assert(registrationIndices.size() == alignmentFileCount);

    std::vector<std::reference_wrapper<const bam_hdr_t>> allHeaders;
    for (unsigned alignmentFileIndex(0); alignmentFileIndex<alignmentFileCount; ++alignmentFileIndex)
    {
        const std::string& alignFile(alignmentFilename[alignmentFileIndex]);
        const unsigned bamIndex(registrationIndices[alignmentFileIndex]);
        const bam_streamer& readStream(streamData.registerBam(alignFile.c_str(), bamIndex));

        allHeaders.push_back(readStream.get_header());

        if (alignmentFileIndex > 0)
        {
            // check that all header chrom details match the first header:
            if (! check_header_compatibility(allHeaders.front(),allHeaders.back()))
            {
                using namespace illumina::common;
                std::ostringstream oss;
                oss << "Input BAM/CRAM files have incompatible headers.\n";
                oss << "\tfile1:\t'" << alignmentFilename.front() << "'\n";
                oss << "\tfile2:\t'" << alignFile << "'\n";
                BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
            }
        }
    }

    return allHeaders;
}



void
getSubRegionsFromBedTrack(
    const std::string& callRegionsBedFilename,
    const std::string& regionChrom,
    const known_pos_range2& regionRange,
    std::vector<known_pos_range2>& subRegionRanges)
{
    static const pos_t maxSubRegionCallGap(5000);

    assert(not callRegionsBedFilename.empty());
    assert(not regionChrom.empty());

    const std::string regionString(getSamtoolsRegionString(regionChrom, regionRange));
    bed_streamer callRegionStream(callRegionsBedFilename.c_str(), regionString.c_str());

    subRegionRanges.clear();
    while (callRegionStream.next())
    {
        bool isStartNewSubRegion(true);
        const bed_record& callRegion(*(callRegionStream.get_record_ptr()));

        // Note that the bed streamer class should, depending on how it was parameterized, either throw an exception or
        // filter invalid records already. Therefore, a simple assertion is sufficient here.
        assert(callRegion.is_valid());

        if (! subRegionRanges.empty())
        {
            {
                // assert that the bed file is sorted (this is redundant with the tabix indexing requirement):
                const pos_t subRegionBegin(subRegionRanges.back().begin_pos());
                assert(std::max(callRegion.begin, regionRange.begin_pos()) >= subRegionBegin);
            }
            const pos_t subRegionEnd(subRegionRanges.back().end_pos());
            const pos_t callGap(callRegion.begin - subRegionEnd);
            isStartNewSubRegion = (callGap > maxSubRegionCallGap);
        }

        const pos_t subRegionEnd = std::min(callRegion.end, regionRange.end_pos());
        if (isStartNewSubRegion)
        {
            const pos_t subRegionBegin = std::max(callRegion.begin, regionRange.begin_pos());
            subRegionRanges.emplace_back(subRegionBegin, subRegionEnd);
        }
        else
        {
            subRegionRanges.back().set_end_pos(std::max(subRegionRanges.back().end_pos(),subRegionEnd));
        }
    }
}



/// This means 'valid' in the sense of what the code can handle right
/// now. Specifically, '=' are not supported.
///
static
bool
is_valid_bam_code(const uint8_t a)
{
    using namespace BAM_BASE;

    switch (a)
    {
    case A:
    case C:
    case G:
    case T:
    case ANY:
        return true;
    default:
        return false;
    }
}



static
bool
is_valid_bam_seq(const bam_seq& bs)
{
    const unsigned rs(bs.size());
    for (unsigned i(0); i<rs; ++i)
    {
        if (! is_valid_bam_code(bs.get_code(i))) return false;
    }
    return true;
}



/// \brief Sanity-check bam record for general validity, and strelka-specific restrictions
///
/// Note that most issues here will trigger an error/exit.
///
/// \return False if read should be filtered
///
static
bool
checkBamRecord(
    const bam_streamer& read_stream,
    const bam_record& read)
{
    const unsigned rs(read.read_size());

    if (rs==0)
    {
        log_os << "ERROR: anomalous read size (<=0) in input alignment record:\n";
        read_stream.report_state(log_os);
        exit(EXIT_FAILURE);
    }

    if (rs > STRELKA_MAX_READ_SIZE)
    {
        log_os << "ERROR: maximum read size (" << STRELKA_MAX_READ_SIZE << ") exceeded in input read alignment record:\n";
        read_stream.report_state(log_os);
        exit(EXIT_FAILURE);
    }

    // check that BAM read sequence contains expected characters:
    const bam_seq bseq(read.get_bam_read());
    if (! is_valid_bam_seq(bseq))
    {
        log_os << "ERROR: unsupported base(s) in read sequence: " << bseq << "\n";
        read_stream.report_state(log_os);
        exit(EXIT_FAILURE);
    }

    // check that BAM qual sequence contains quality values we can handle:
    {
        const uint8_t* qual(read.qual());
        for (unsigned i(0); i<rs; ++i)
        {
            try
            {
                // I *think* this is here because a QUAL value of "*" (in SAM), will be translated to
                // a basecall quality of 255 to indicate that the basecall quality is unknown.
                //
                // There was also an older RNAseq work-flow corner case where reads of length one with '*' QUAL
                // fields were used as part of the RNA mapper's output scheme.
                //
                if (qual[i]==255)
                {
                    return false;
                }
                qphred_cache::qscore_check(qual[i],"basecall quality");
            }
            catch (...)
            {
                log_os << "\nException for basecall quality score " << static_cast<int>(qual[i]) << " at read position " << (i+1) << "\n";
                read_stream.report_state(log_os);
                throw;
            }
        }
    }

    return true;
}



static
bool
is_usable_read_mapping(const starling_base_options& opt,
                       const bam_record& read,
                       const bool is_tier2 =false)
{
    int current_min_mapping_quality(opt.minMappingErrorPhredProb);
    bool is_include_singleton(opt.includeSingletonReads);
    bool is_include_anomalous(opt.includeAnomalousReads);
    if (is_tier2)
    {
        current_min_mapping_quality=opt.tier2.minMappingErrorPhredProb;
        if (opt.tier2.includeSingletonReads)
        {
            is_include_singleton=true;
        }
        if (opt.tier2.includeAnomalousReads)
        {
            is_include_anomalous=true;
        }
    }

    if (read.is_paired())
    {
        const bool is_singleton(read.is_mate_unmapped());
        const bool is_anomalous((! is_singleton) && (! read.is_proper_pair()));
        if        ((! is_include_singleton) && is_singleton)
        {
            return false;
        }
        else if ((! is_include_anomalous) && is_anomalous)
        {
            return false;
        }
    }

    {
        const int mapq(static_cast<int>(read.map_qual()));
        if (mapq < current_min_mapping_quality)
        {
            return false;
        }
    }
    return true;
}



static
MAPLEVEL::index_t
get_map_level(const starling_base_options& opt,
              const bam_record& read)
{
    using namespace MAPLEVEL;

    if (read.is_unmapped()) return UNMAPPED;

    if (is_usable_read_mapping(opt,read)) return TIER1_MAPPED;
    if (opt.useTier2Evidence)
    {
        if (is_usable_read_mapping(opt,read,true)) return TIER2_MAPPED;
    }
    return SUB_MAPPED;
}



/// check if read is overdepth, but only if max input depth option is set
static
bool
isAlignmentAboveMaxDepth(
    const starling_base_options& opt,
    const starling_pos_processor_base& posProcessor,
    const unsigned sample_no,
    const alignment& al)
{
    using namespace ALIGNPATH;

    if (! opt.is_max_input_depth) return false;

    pos_t ref_head_pos(al.pos);

    for (const auto& ps : al.path)
    {
        if (is_segment_align_match(ps.type))
        {
            if (posProcessor.is_estimated_depth_range_ge_than(ref_head_pos,
                                                              ref_head_pos+static_cast<pos_t>(ps.length),
                                                              opt.max_input_depth,
                                                              sample_no))
            {
                return true;
            }
        }

        if (is_segment_type_ref_length(ps.type)) ref_head_pos += ps.length;
    }
    return false;
}



void
processInputReadAlignment(
    const starling_base_options& opt,
    const reference_contig_segment& ref,
    const bam_streamer& read_stream,
    const bam_record& read,
    const pos_t base_pos,
    starling_read_counts& readCounts,
    starling_pos_processor_base& posProcessor,
    const unsigned sampleIndex)
{
    // These read filters are always on, because we can't do anything sensible
    // with these cases (except supplement). This filtration is shared with the chromosome
    // depth estimation routine, so the perceived depth of the two routines
    // will match for the purpose of high/low depth filtration.
    //
    const READ_FILTER_TYPE::index_t filterIndex(starling_read_filter_shared(read));
    if (filterIndex != READ_FILTER_TYPE::NONE)
    {
        using namespace READ_FILTER_TYPE;
        if (filterIndex == PRIMARY) readCounts.primary_filter++;
        if (filterIndex == DUPLICATE) readCounts.duplicate++;
        if (filterIndex == UNMAPPED) readCounts.unmapped++;
        if (filterIndex == SECONDARY) readCounts.secondary++;
        if (filterIndex == SUPPLEMENTARY) readCounts.supplement++;
        return;
    }

    MAPLEVEL::index_t maplev(get_map_level(opt,read));

    const bool isKeepRecord(checkBamRecord(read_stream, read));
    if (! isKeepRecord)
    {
        return;
    }

    // Now that 'alignment-free' filtration is finished, sanity-check/filter/normalize the remaining mapped reads
    // based on the read's alignment.
    //
    // Note there is logic below assuming that unmapped reads need to be handled even though these
    // are filtered out above. This is to retain the ability to experiment with local assembly of
    // unmapped reads with mapped mate reads.
    //

    // if we're seeing an unmapped read, it's still expected to have a position and chrom assignment:
    alignment readAlignment;
    readAlignment.pos=base_pos;

    if (maplev != MAPLEVEL::UNMAPPED)
    {
        // if mapped, sanity check alignment:
        readAlignment.is_fwd_strand=read.is_fwd_strand();
        bam_cigar_to_apath(read.raw_cigar(),read.n_cigar(),readAlignment.path);

        ALIGNPATH::apath_cleaner(readAlignment.path);

        if (read.read_size() != ALIGNPATH::apath_read_length(readAlignment.path))
        {
            using namespace illumina::common;

            const unsigned rs(read.read_size());
            const unsigned as(ALIGNPATH::apath_read_length(readAlignment.path));

            std::ostringstream oss;
            oss << "Read length implied by mapped alignment (" << as << ") does not match read length ("
                << rs << ") in alignment record:\n";
            read_stream.report_state(oss);
            BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
        }

        // filter out reads with no match segments:
        if (ALIGNPATH::is_apath_floating(readAlignment.path))
        {
            readCounts.floating++;
            return;
        }

        if (isAlignmentAboveMaxDepth(opt, posProcessor, sampleIndex, readAlignment))
        {
            readCounts.max_depth++;
            return;
        }

        // normalize/left-shift the input alignment
        const rc_segment_bam_seq refBamSeq(ref);
        const bam_seq readBamSeq(read.get_bam_read());
        normalizeAlignment(refBamSeq, readBamSeq, readAlignment);
    }


    try
    {
        const char* chrom_name(read_stream.target_id_to_name(read.target_id()));
        posProcessor.insert_read(read,readAlignment,chrom_name,maplev,sampleIndex);
    }
    catch (...)
    {
        log_os << "\nException caught while inserting read alignment in posProcessor. Genomic read alignment record:\n";
        read_stream.report_state(log_os);
        throw;
    }
}



/// return common prefix and suffix length of two strings, where
/// suffix takes priority (this is important for left-shifting)
///
static
std::pair<unsigned,unsigned>
common_xfix_length(
    const std::string& s1,
    const std::string& s2)
{
    const unsigned s1s(s1.size());
    const unsigned s2s(s2.size());

    unsigned suffix(0);
    unsigned nsearch(std::min(s1s,s2s));
    for (; suffix<nsearch; ++suffix)
    {
        if (s1[s1s-1-suffix] != s2[s2s-1-suffix]) break;
    }

    unsigned prefix(0);
    nsearch -= suffix;
    for (; prefix<nsearch; ++prefix)
    {
        if (s1[prefix] != s2[prefix]) break;
    }

    return std::make_pair(prefix,suffix);
}



bool
convert_vcfrecord_to_indel_allele(
    const unsigned max_indel_size,
    const vcf_record& vcf_indel,
    const unsigned altIndex,
    IndelObservation& obs)
{
    assert(vcf_indel.is_indel());
    assert(altIndex<vcf_indel.alt.size());

    const unsigned rs(vcf_indel.ref.size());
    const auto& alt(vcf_indel.alt[altIndex]);
    const unsigned as(alt.size());
    const std::pair<unsigned,unsigned> xfix(common_xfix_length(vcf_indel.ref,alt));
    const unsigned nfix(xfix.first+xfix.second);
    assert(nfix<=std::min(rs,as));
    const int insert_length(as-nfix);
    const int delete_length(rs-nfix);

    if ((insert_length > static_cast<int>(max_indel_size)) ||
        (delete_length > static_cast<int>(max_indel_size))) return false;

    if ((insert_length==0) and (delete_length==0))
    {
        using namespace illumina::common;
        std::ostringstream oss;
        oss << "Can't parse vcf indel: '" << vcf_indel << "'";
        BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
    }

    // starling indel pos is at the first changed base but zero-indexed:
    obs.key.pos = (vcf_indel.pos+xfix.first-1);
    obs.key.type = INDEL::INDEL;
    obs.key.deletionLength = delete_length;
    obs.key.insertSequence = std::string(alt.begin()+xfix.first,alt.end()-xfix.second);
    return true;
}



void
process_candidate_indel(
    const unsigned max_indel_size,
    const vcf_record& vcf_indel,
    starling_pos_processor_base& posProcessor,
    const unsigned sampleIndex,
    const bool is_forced_output)
{
    assert (vcf_indel.is_indel());

    const unsigned altCount(vcf_indel.alt.size());
    for (unsigned altIndex(0); altIndex<altCount; ++altIndex)
    {
        IndelObservation obs;
        const bool isAlleleConverted =
            convert_vcfrecord_to_indel_allele(max_indel_size,vcf_indel,altIndex,obs);
        if (! isAlleleConverted) continue;

        obs.data.is_external_candidate = true;
        obs.data.is_forced_output = is_forced_output;

        posProcessor.insert_indel(obs,sampleIndex);
    }
}
