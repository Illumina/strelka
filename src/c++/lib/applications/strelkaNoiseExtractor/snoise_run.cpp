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

#include "snoise_run.hh"
#include "snoise_pos_processor.hh"
#include "snoise_streams.hh"

#include "blt_util/log.hh"
#include "common/Exceptions.hh"
#include "htsapi/bam_header_info.hh"
#include "starling_common/starling_base_shared.hh"
#include "starling_common/HtsMergeStreamerUtil.hh"
#include "starling_common/starling_ref_seq.hh"
#include "starling_common/starling_pos_processor_util.hh"



namespace INPUT_TYPE
{
enum index_t
{
    CANDIDATE_INDELS,
    FORCED_GT_VARIANTS,
};
}



void
snoise_run(
    const prog_info& pinfo,
    const snoise_options& opt)
{
    starling_read_counts brc;
    reference_contig_segment ref;

    ////////////////////////////////////////
    // setup streamData:
    //
    HtsMergeStreamer streamData;

    // additional data structures required in the region loop below, which are filled in as a side effect of
    // streamData initialization:
    std::vector<std::reference_wrapper<const bam_hdr_t>> bamHeaders;
    {
        std::vector<unsigned> registrationIndices(opt.alignFileOpt.alignmentFilename.size(), 0);
        bamHeaders = registerAlignments(opt.alignFileOpt.alignmentFilename, registrationIndices, streamData);

        assert(not bamHeaders.empty());
        const bam_hdr_t& referenceHeader(bamHeaders.front());

        registerVcfList(opt.input_candidate_indel_vcf, INPUT_TYPE::CANDIDATE_INDELS, referenceHeader, streamData);
        registerVcfList(opt.force_output_vcf, INPUT_TYPE::FORCED_GT_VARIANTS, referenceHeader, streamData);
    }

    const bam_hdr_t& referenceHeader(bamHeaders.front());
    const bam_header_info referenceHeaderInfo(referenceHeader);

    const unsigned regionCount(opt.regions.size());
    for (unsigned regionIndex(0); regionIndex<regionCount; ++regionIndex)
    {
        const std::string& region(opt.regions[regionIndex]);
        AnalysisRegionInfo rinfo;
        getStrelkaAnalysisRegions(region, opt.max_indel_size, rinfo);

        // check that target region chrom exists in bam headers:
        if (not referenceHeaderInfo.chrom_to_index.count(rinfo.regionChrom))
        {
            using namespace illumina::common;
            std::ostringstream oss;
            oss << "ERROR: region contig name: '" << rinfo.regionChrom << "' is not found in the header of BAM/CRAM file: '" << opt.alignFileOpt.alignmentFilename.front() << "'\n";
            BOOST_THROW_EXCEPTION(LogicException(oss.str()));
        }

        streamData.resetRegion(rinfo.streamerRegion.c_str());
        setRefSegment(opt, rinfo.regionChrom, rinfo.refRegionRange, ref);

    const starling_base_deriv_options dopt(opt, ref);
    const unsigned sampleCount(1);
    snoise_streams streams(opt, pinfo, referenceHeader, sampleCount);

    snoise_pos_processor sppr(opt, dopt, ref, streams);

        sppr.resetRegion(rinfo.regionChrom, rinfo.regionRange);

    while (streamData.next())
    {
        const pos_t currentPos(streamData.getCurrentPos());
        const HTS_TYPE::index_t currentHtsType(streamData.getCurrentType());
        const unsigned currentIndex(streamData.getCurrentIndex());

        if (currentPos >= rinfo.streamerRegionRange.end_pos()) break;

        // wind sppr forward to position behind buffer head:
        sppr.set_head_pos(currentPos - 1);

        if (HTS_TYPE::BAM == currentHtsType)
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
                                      brc, sppr);
        }
        else if (HTS_TYPE::VCF == currentHtsType)
        {
            const vcf_record& vcfRecord(streamData.getCurrentVcf());
            if (INPUT_TYPE::CANDIDATE_INDELS == currentIndex)     // process candidate indels input from vcf file(s)
            {
                if (vcfRecord.is_indel())
                {
                    process_candidate_indel(opt.max_indel_size, vcfRecord, sppr);
                }
            }
            else if (INPUT_TYPE::FORCED_GT_VARIANTS ==
                     currentIndex)     // process forced genotype tests from vcf file(s)
            {
                if (vcfRecord.is_indel())
                {
                    static const unsigned sample_no(0);
                    static const bool is_forced_output(true);
                    process_candidate_indel(opt.max_indel_size, vcfRecord, sppr, sample_no, is_forced_output);
                }
                else if (vcfRecord.is_snv())
                {
                    sppr.insert_forced_output_pos(vcfRecord.pos - 1);
                }
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
}
