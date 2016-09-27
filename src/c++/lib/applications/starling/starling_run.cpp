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

#include <blt_util/id_map.hh>
#include "starling_run.hh"
#include "starling_pos_processor.hh"
#include "starling_streams.hh"

#include "appstats/RunStatsManager.hh"
#include "blt_util/log.hh"
#include "common/Exceptions.hh"
#include "starling_common/HtsMergeStreamerUtil.hh"
#include "starling_common/ploidy_util.hh"
#include "starling_common/starling_ref_seq.hh"
#include "starling_common/starling_pos_processor_util.hh"



namespace INPUT_TYPE
{
enum index_t
{
    CANDIDATE_INDELS,
    FORCED_GT_VARIANTS,
    PLOIDY_REGION,
    NOCOMPRESS_REGION
};
}



/// (1) validate sample count/sample name match to input alignment files
/// (2) construct a sample index translation map
///
/// \param[out] sampleIndexToPloidyVcfSampleIndex translate from the primary sample index order
///                 (derived from input alignment file order) to the VCF sample index (index in VCF sample columns)
///
/// TODO better place to move this fn?
static
void
mapVcfSampleIndices(
    const vcf_streamer& vcfStream,
    const std::vector<std::string>& sampleNames,
    std::vector<unsigned> sampleIndexToPloidyVcfSampleIndex)
{
    using namespace illumina::common;

    id_set<std::string> vcfSampleSet;
    const unsigned vcfSampleCount(vcfStream.getSampleCount());
    for (unsigned vcfSampleIndex(0); vcfSampleIndex < vcfSampleCount; ++vcfSampleIndex)
    {
        const std::string vcfSampleName(vcfStream.getSampleName(vcfSampleIndex));

        if (vcfSampleSet.test_key(vcfSampleName))
        {
            std::ostringstream oss;
            oss << "ERROR: repeated entry for sample name '" << vcfSampleName << "' in ploidy VCF file: '" << vcfStream.name() << "'\n";
            BOOST_THROW_EXCEPTION(LogicException(oss.str()));
        }

        vcfSampleSet.insert_key(vcfSampleName);
    }

    for (const auto& sampleName : sampleNames)
    {
        const auto maybeId(vcfSampleSet.get_optional_id(sampleName));
        if (not maybeId)
        {
            std::ostringstream oss;
            oss << "ERROR: no entry for sample name '" << sampleName << "' in ploidy VCF file: '" << vcfStream.name() << "'\n";
            BOOST_THROW_EXCEPTION(LogicException(oss.str()));
        }
        else
        {
            sampleIndexToPloidyVcfSampleIndex.push_back(*maybeId);
        }
    }
}



void
starling_run(
    const prog_info& pinfo,
    const starling_options& opt)
{
    opt.validate();

    RunStatsManager segmentStatMan(opt.segmentStatsFilename);

    reference_contig_segment ref;
    get_starling_ref_seq(opt,ref);

    const starling_deriv_options dopt(opt,ref);
    const pos_range& rlimit(dopt.report_range_limit);

    const std::string bam_region(get_starling_bam_region_string(opt,dopt));

    HtsMergeStreamer streamData(bam_region.c_str());

    const unsigned sampleCount(opt.alignFileOpt.alignmentFilename.size());

    std::vector<unsigned> registrationIndices;
    for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
    {
        registrationIndices.push_back(sampleIndex);
    }
    const auto allHeaders(registerAlignments(opt, opt.alignFileOpt, registrationIndices, streamData));
    starling_streams client_io(opt, pinfo, allHeaders, sampleCount);

    starling_pos_processor sppr(opt,dopt,ref,client_io);
    starling_read_counts brc;

    assert(not allHeaders.empty());
    const bam_hdr_t& referenceHeader(allHeaders.front());

    static const bool noRequireNormalized(false);
    registerVcfList(opt.input_candidate_indel_vcf, INPUT_TYPE::CANDIDATE_INDELS, referenceHeader, streamData, noRequireNormalized);
    registerVcfList(opt.force_output_vcf, INPUT_TYPE::FORCED_GT_VARIANTS, referenceHeader, streamData);

    std::vector<unsigned> sampleIndexToPloidyVcfSampleIndex;
    if (! opt.ploidy_region_vcf.empty())
    {
        const vcf_streamer& vcfStream(streamData.registerVcf(opt.ploidy_region_vcf.c_str(), INPUT_TYPE::PLOIDY_REGION));
        vcfStream.validateBamHeaderChromSync(referenceHeader);

        const std::vector<std::string>& sampleNames();
        mapVcfSampleIndices(vcfStream, client_io.getSampleNames(), sampleIndexToPloidyVcfSampleIndex);
    }

    if (! opt.gvcf.nocompress_region_bedfile.empty())
    {
        streamData.registerBed(opt.gvcf.nocompress_region_bedfile.c_str(), INPUT_TYPE::NOCOMPRESS_REGION);
    }

    while (streamData.next())
    {
        const pos_t currentPos(streamData.getCurrentPos());
        const HTS_TYPE::index_t currentHtsType(streamData.getCurrentType());
        const unsigned currentIndex(streamData.getCurrentIndex());

        // Process finishes at the the end of rlimit range. Note that
        // some additional padding is allowed for off-range indels
        // which might influence results within rlimit:
        //
        if (rlimit.is_end_pos && (currentPos >= (rlimit.end_pos+static_cast<pos_t>(opt.max_indel_size)))) break;

        // wind sppr forward to position behind buffer head:
        sppr.set_head_pos(currentPos-1);

        if       (HTS_TYPE::BAM == currentHtsType)
        {
            // Remove the filter below because it's not valid for
            // RNA-Seq case, reads should be selected for the report
            // range by the bam reading functions
            //
            // /// get potential bounds of the read based only on current_pos:
            // const known_pos_range any_read_bounds(current_pos-max_indel_size,current_pos+MAX_READ_SIZE+max_indel_size);
            // if( sppr.is_range_outside_report_influence_zone(any_read_bounds) ) continue;

            // Approximate begin range filter: (removed for RNA-Seq)
            //if((current_pos+MAX_READ_SIZE+max_indel_size) <= rlimit.begin_pos) continue;

            processInputReadAlignment(opt, ref, streamData.getCurrentBamStreamer(),
                                      streamData.getCurrentBam(), currentPos,
                                      rlimit.begin_pos, brc, sppr, currentIndex);
        }
        else if (HTS_TYPE::VCF == currentHtsType)
        {
            const vcf_record& vcfRecord(streamData.getCurrentVcf());
            if     (INPUT_TYPE::CANDIDATE_INDELS == currentIndex)     // process candidate indels input from vcf file(s)
            {
                if (vcfRecord.is_indel())
                {
                    process_candidate_indel(opt.max_indel_size, vcfRecord, sppr);
                }
            }
            else if (INPUT_TYPE::FORCED_GT_VARIANTS == currentIndex)     // process forced genotype tests from vcf file(s)
            {
                if (vcfRecord.is_indel())
                {
                    static const unsigned sample_no(0);
                    static const bool is_forced_output(true);
                    process_candidate_indel(opt.max_indel_size, vcfRecord,sppr,sample_no,is_forced_output);
                }
                else if (vcfRecord.is_snv())
                {
                    sppr.insert_forced_output_pos(vcfRecord.pos-1);
                }
            }
            else
            {
                assert(false && "Unexpected hts index");
            }
        }
        else if (HTS_TYPE::BED == currentHtsType)
        {
            const bed_record& bedRecord(streamData.getCurrentBed());
            if     (INPUT_TYPE::PLOIDY_REGION == currentIndex)
            {
                known_pos_range2 ploidyRange(bedRecord.begin,bedRecord.end);
                const unsigned ploidy(parsePloidyFromBedStrict(bedRecord.line));
                if ((ploidy == 0) || (ploidy == 1))
                {
                    const bool retval(sppr.insert_ploidy_region(ploidyRange,ploidy));
                    if (! retval)
                    {
                        using namespace illumina::common;

                        std::ostringstream oss;
                        oss << "ERROR: ploidy bedfile record conflicts with prior record. Bedfile line: '" << bedRecord.line << "'\n";
                        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
                    }
                }
            }
            else if (INPUT_TYPE::NOCOMPRESS_REGION == currentIndex)
            {
                known_pos_range2 range(bedRecord.begin,bedRecord.end);
                sppr.insert_nocompress_region(range);
            }
            else
            {
                assert(false && "Unexpected hts index");
            }
        }
        else
        {
            log_os << "ERROR: invalid input condition.\n";
            exit(EXIT_FAILURE);
        }
    }

    sppr.reset();
}

