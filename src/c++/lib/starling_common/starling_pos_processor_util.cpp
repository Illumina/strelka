// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
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

///
/// \author Chris Saunders
///
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the beginning of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///


#include "candidate_alignment.hh"
#include "normalizeAlignment.hh"
#include "starling_pos_processor_indel_util.hh"
#include "starling_read_filter_shared.hh"
#include "starling_read_util.hh"

#include "blt_util/align_path.hh"
#include "blt_util/log.hh"
#include "blt_util/qscore.hh"
#include "common/Exceptions.hh"
#include "htsapi/align_path_bam_util.hh"
#include "htsapi/bam_seq.hh"
#include "starling_common/starling_pos_processor_util.hh"

#include <cassert>

#include <sstream>



std::string
get_starling_bam_region_string(const starling_base_options& opt,
                               const starling_base_deriv_options& dopt)
{
    const int zsize(opt.max_indel_size);
    const pos_t begin_pos(std::max(0,dopt.report_range.begin_pos-zsize));
    const pos_t end_pos(dopt.report_range.end_pos+zsize);

    std::ostringstream bam_region_oss;
    bam_region_oss << opt.bam_seq_name << ':' << begin_pos+1 << '-' << end_pos;
    return bam_region_oss.str();
}





// This means 'valid' in the sense of what the code can handle right
// now. Specifically, '=' are not supported.
//
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

static
bool
check_bam_record(const bam_streamer& read_stream,
                 const bam_record& read)
{
    const unsigned rs(read.read_size());

    if (rs==0)
    {
        log_os << "ERROR: anomalous read size (<=0) in input alignment record:\n";
        read_stream.report_state(log_os);
        exit(EXIT_FAILURE);
    }

    if (rs > STARLING_MAX_READ_SIZE)
    {
        log_os << "ERROR: maximum read size (" << STARLING_MAX_READ_SIZE << ") exceeded in input read alignment record:\n";
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
                //for RNAseq work-flow corner case. Reads of length one with '*' q-score
//                log_os << "qscore " << static_cast<int>(qual[i])<< "\n";
                if (qual[i]==255)
                {
//                    log_os << "\nException for basecall quality score " << static_cast<int>(qual[i])<< "\n";
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
    int current_min_mapping_quality(opt.min_mapping_quality);
    bool is_filter_unanchored(opt.is_filter_unanchored);
    bool is_include_singleton(opt.is_include_singleton);
    bool is_include_anomalous(opt.is_include_anomalous);
    if (is_tier2)
    {
        if (opt.tier2.is_tier2_min_mapping_quality)
        {
            current_min_mapping_quality=opt.tier2.tier2_min_mapping_quality;
        }
        if (opt.tier2.is_tier2_no_filter_unanchored)
        {
            is_filter_unanchored=false;
        }
        if (opt.tier2.is_tier2_include_singleton)
        {
            is_include_singleton=true;
        }
        if (opt.tier2.is_tier2_include_anomalous)
        {
            is_include_anomalous=true;
        }
    }

    if (is_filter_unanchored &&
        read.is_unanchored())
    {
        return false;
    }

    if (read.is_paired())
    {
        const bool is_singleton(read.is_mate_unmapped());
        const bool is_anomalous((! is_singleton) && (! read.is_proper_pair()));
        if        ((! is_include_singleton) && is_singleton)
        {
            return false; // singleton
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
    if (opt.tier2.is_tier2())
    {
        if (is_usable_read_mapping(opt,read,true)) return TIER2_MAPPED;
    }
    return SUB_MAPPED;
}



static
bool
is_al_overdepth(const starling_base_options& opt,
                const starling_pos_processor_base& sppr,
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
            if (sppr.is_estimated_depth_range_ge_than(ref_head_pos,
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



// handles genomic read alignments -- reads are parsed, their indels
// are extracted and buffered, and the reads themselves are buffered
//
void
process_genomic_read(
    const starling_base_options& opt,
    const reference_contig_segment& ref,
    const bam_streamer& read_stream,
    const bam_record& read,
    const pos_t base_pos,
    const pos_t /*report_begin_pos*/,
    starling_read_counts& brc,
    starling_pos_processor_base& sppr,
    const unsigned sample_no)
{
    // Read filters which are *always* on, because Strelka
    // can't do anything sensible with these reads. This filtration
    // is shared with the chrom depth estimation routine.
    //
    const READ_FILTER_TYPE::index_t filterId(starling_read_filter_shared(read));
    if (filterId != READ_FILTER_TYPE::NONE)
    {
        using namespace READ_FILTER_TYPE;
        if (filterId == PRIMARY) brc.primary_filter++;
        if (filterId == DUPLICATE) brc.duplicate++;
        if (filterId == UNMAPPED) brc.unmapped++;
        if (filterId == SECONDARY) brc.secondary++;
        if (filterId == SUPPLEMENT) brc.supplement++;
        return;
    }

    MAPLEVEL::index_t maplev(get_map_level(opt,read));

    // RNAseq modification, qscore 255 used as flag
    // do not throw exception but rather ignore read
    // logic implemented in check_bam_record
    const bool allGood  = check_bam_record(read_stream,read);
    if (!allGood)
    {
//        log_os << "Skipping read " << read << "\n";
        return;
    }

    // secondary range filter check:
    //   (removed as part of RNA-Seq modifications)
    //
    //if(base_pos+static_cast<pos_t>(read.read_size())+MAX_INDEL_SIZE <= report_begin_pos) return;



    // extract indels and add to indel map:
    //
    // TODO: settle how to handle reads which fail mapping score but could be realigned
    //
    // include submapped reads for realignment only if we are writing out the realigned cases
    // ** feature disabled for now -- unmapped reads are ignored unless they come from grouper
    //
    {
        //    if(is_usable_mapping or opt.is_realign_submapped_reads){ // or opt.is_realigned_read_file){

        // if we're seeing an unmapped read, it's still expected to
        // have a position and chrom assignment
        //
        alignment al;
        al.pos=base_pos;

        // if mapped, sanity check alignment:
        if (maplev != MAPLEVEL::UNMAPPED)
        {
            al.is_fwd_strand=read.is_fwd_strand();
            bam_cigar_to_apath(read.raw_cigar(),read.n_cigar(),al.path);

            ALIGNPATH::apath_cleaner(al.path);

            if (read.read_size() != ALIGNPATH::apath_read_length(al.path))
            {
                const unsigned rs(read.read_size());
                const unsigned as(ALIGNPATH::apath_read_length(al.path));
                log_os << "ERROR: Read length implied by mapped alignment (" << as << ") does not match read length (" << rs << ") in alignment record:\n";
                read_stream.report_state(log_os);
                exit(EXIT_FAILURE);
            }

            // filter out reads with no match segments:
            if (ALIGNPATH::is_apath_floating(al.path))
            {
                brc.floating++;
                return;
            }

            if (is_al_overdepth(opt,sppr,sample_no,al))
            {
                brc.max_depth++;
                return;
            }

            // normalize/left-shift the input alignment before inserting into the read buffer below:
            const rc_segment_bam_seq refBamSeq(ref);
            const bam_seq readBamSeq(read.get_bam_read());
            normalizeAlignment(refBamSeq,readBamSeq,al);
        }


        try
        {
            const char* chrom_name(read_stream.target_id_to_name(read.target_id()));
            sppr.insert_read(read,al,chrom_name,maplev,sample_no);
        }
        catch (...)
        {
            log_os << "\nException caught while inserting read alignment in sppr. Genomic read alignment record:\n";
            read_stream.report_state(log_os);
            throw;
        }
    }
}



// return common prefix and suffix length of two strings, where
// suffix takes priority (this is important for left-shifting)
//
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
        oss << "ERROR: Can't parse vcf indel: '" << vcf_indel << "'\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    // starling indel pos is at the first changed base but zero-indexed:
    obs.key.pos = (vcf_indel.pos+xfix.first-1);
    obs.key.type = INDEL::INDEL;
    obs.key.deletionLength = delete_length;
    obs.key.insertSequence = std::string(alt.begin()+xfix.first,alt.end()-xfix.second);
    return true;
}



// handles candidate indel input from a vcf record
//
void
process_candidate_indel(
    const unsigned max_indel_size,
    const vcf_record& vcf_indel,
    starling_pos_processor_base& sppr,
    const unsigned sample_no,
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

        sppr.insert_indel(obs,sample_no);
    }
}
