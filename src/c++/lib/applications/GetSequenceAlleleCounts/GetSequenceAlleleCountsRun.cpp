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

#include "GetSequenceAlleleCountsRun.hh"
#include "SequenceAlleleCountsPosProcessor.hh"

#include "appstats/RunStatsManager.hh"
#include "blt_util/log.hh"
#include "common/Exceptions.hh"
#include "htsapi/align_path_bam_util.hh"
#include "htsapi/bam_header_info.hh"
#include "htsapi/vcf_record_util.hh"
#include "starling_common/HtsMergeStreamerUtil.hh"
#include "starling_common/ploidy_util.hh"
#include "starling_common/starling_ref_seq.hh"
#include "starling_common/starling_pos_processor_util.hh"
#include "strelka_common/StrelkaSampleSetSummary.hh"



namespace INPUT_TYPE
{
enum index_t
{
    CANDIDATE_INDELS,
    FORCED_GT_VARIANTS,
    KNOWN_VARIANTS,
    EXCLUDE_REGION,
};
}


/// \brief Adds truth variant input from a vcf into the error counting process
///
static
void
processTrueIndelVariantRecord(
    const unsigned max_indel_size,
    const vcf_record& vcfr,
    SequenceAlleleCountsPosProcessor& posProcessor)
{
    /// TODO: the actual truth logic
    ///
    /// to make the problem incremental/easier consider:
    /// 1) only handle indels at first
    /// 2) ignore any record with multiple alts.
    ///

    // N.B. ignoring SNPs and indels with multiple alts for now
    if (vcfr.is_indel() && vcfr.alt.size() == 1)
    {
        const unsigned altCount(vcfr.alt.size());

        for (unsigned altIndex(0); altIndex<altCount; ++altIndex)
        {
            IndelObservation obs;
            const bool isAlleleConverted =
                convert_vcfrecord_to_indel_allele(max_indel_size,vcfr,altIndex,obs);
            if (! isAlleleConverted) continue;

            obs.data.is_external_candidate = true;
            // obs.data.is_forced_output = is_forced_output;
            // not setting is_forced_output since it's not clear what that would mean in
            // the pattern analyzer

            const unsigned sample_no(0); //currently, the pattern analyzer only operates on a single sample

            posProcessor.insert_indel(obs,sample_no);
            posProcessor.addKnownVariant(vcfr);
        }
    }
}



void
getSequenceAlleleCountsRun(
    const prog_info& pinfo,
    const SequenceAlleleCountsOptions& opt)
{
    // ensure that this object is created first to improve accuracy of runtime benchmarking
    RunStatsManager statsManager(opt.segmentStatsFilename);

    opt.validate();

    const SequenceAlleleCountsDerivOptions dopt(opt);
    starling_read_counts readCounts;
    reference_contig_segment ref;

    ////////////////////////////////////////
    // setup streamData:
    //
    HtsMergeStreamer streamData(opt.referenceFilename);

    // additional data structures required in the region loop below, which are filled in as a side effect of
    // streamData initialization:
    std::vector<std::reference_wrapper<const bam_hdr_t>> bamHeaders;
    {
        if (opt.alignFileOpt.alignmentFilenames.size() > 1)
        {
            using namespace illumina::common;
            std::ostringstream oss;
            oss << "WARNING: Multiple bam file inputs. Will be treated as single sample (with sampleName: " << opt.alignFileOpt.alignmentFilenames[0] << ")\n";
        }
        std::vector<unsigned> registrationIndices(opt.alignFileOpt.alignmentFilenames.size(), 0);
        bamHeaders = registerAlignments(opt.alignFileOpt.alignmentFilenames, registrationIndices, streamData);

        assert(! bamHeaders.empty());
        const bam_hdr_t& referenceHeader(bamHeaders.front());

        static const bool noRequireNormalized(false);
        registerVcfList(opt.input_candidate_indel_vcf, INPUT_TYPE::CANDIDATE_INDELS, referenceHeader,
                        streamData, noRequireNormalized);
        registerVcfList(opt.force_output_vcf, INPUT_TYPE::FORCED_GT_VARIANTS, referenceHeader, streamData);

        if (!opt.knownVariantsFile.empty())
        {
            const vcf_streamer& vcfStream(
                streamData.registerVcf(opt.knownVariantsFile.c_str(), INPUT_TYPE::KNOWN_VARIANTS));
            vcfStream.validateBamHeaderChromSync(referenceHeader);
        }

        for (const std::string& excludeRegionFilename : opt.excludedRegionsFileList)
        {
            streamData.registerBed(excludeRegionFilename.c_str(), INPUT_TYPE::EXCLUDE_REGION);
        }
    }

    const bam_hdr_t& referenceHeader(bamHeaders.front());
    const bam_header_info referenceHeaderInfo(referenceHeader);

    SequenceAlleleCountsStreams fileStreams(opt, pinfo, referenceHeader);
    SequenceAlleleCountsPosProcessor posProcessor(opt, dopt, ref, fileStreams, statsManager);

    // parse and sanity check regions
    assert ((! opt.isHaplotypingEnabled) && "Region border size must be updated if haplotyping is enabled");
    const unsigned supplementalRegionBorderSize(opt.maxIndelSize);

    const auto& referenceAlignmentFilename(opt.alignFileOpt.alignmentFilenames.front());
    std::vector<AnalysisRegionInfo> regionInfoList;
    getStrelkaAnalysisRegions(opt, referenceAlignmentFilename, referenceHeaderInfo, supplementalRegionBorderSize, regionInfoList);

    for (const auto& regionInfo : regionInfoList)
    {
        posProcessor.resetRegion(regionInfo.regionChrom, regionInfo.regionRange);
        streamData.resetRegion(regionInfo.streamerRegion.c_str());
        setRefSegment(opt, regionInfo.regionChrom, regionInfo.refRegionRange, ref);

        while (streamData.next())
        {
            const pos_t currentPos(streamData.getCurrentPos());
            const HTS_TYPE::index_t currentHtsType(streamData.getCurrentType());
            const unsigned currentIndex(streamData.getCurrentIndex());

            if (currentPos >= regionInfo.streamerRegionRange.end_pos()) break;

            // wind posProcessor forward to position behind buffer head:
            posProcessor.set_head_pos(currentPos - 1);

            if (HTS_TYPE::BAM == currentHtsType)
            {
                // Remove the filter below because it's not valid for
                // RNA-Seq case, reads should be selected for the report
                // range by the bam reading functions
                //
                // /// get potential bounds of the read based only on current_pos:
                // const known_pos_range any_read_bounds(current_pos-maxIndelSize,current_pos+MAX_READ_SIZE+maxIndelSize);
                // if( posProcessor.is_range_outside_report_influence_zone(any_read_bounds) ) continue;

                // Approximate begin range filter: (removed for RNA-Seq)
                //if((current_pos+MAX_READ_SIZE+maxIndelSize) <= rlimit.begin_pos) continue;

                const bam_record& read(streamData.getCurrentBam());

                // special read noise filter used in error counting only -- this isn't the ideal place for this logic:
                {
                    using namespace ALIGNPATH;
                    path_t apath;
                    bam_cigar_to_apath(read.raw_cigar(), read.n_cigar(), apath);

                    if (apath_indel_count(apath) > 2) continue;
                }

                processInputReadAlignment(opt, ref, streamData.getCurrentBamStreamer(),
                                          read, currentPos, readCounts, posProcessor);
            }
            else if (HTS_TYPE::VCF == currentHtsType)
            {
                assertExpectedVcfReference(ref, streamData.getCurrentVcfStreamer());
                const vcf_record& vcfRecord(streamData.getCurrentVcf());
                if (INPUT_TYPE::CANDIDATE_INDELS == currentIndex)     // process candidate indels input from vcf file(s)
                {
                    if (vcfRecord.is_indel())
                    {
                        process_candidate_indel(opt.maxIndelSize, vcfRecord, posProcessor);
                    }
                }
                else if (INPUT_TYPE::FORCED_GT_VARIANTS ==
                         currentIndex)     // process forced genotype tests from vcf file(s)
                {
                    if (vcfRecord.is_indel())
                    {
                        static const unsigned sample_no(0);
                        static const bool is_forced_output(true);
                        process_candidate_indel(opt.maxIndelSize, vcfRecord, posProcessor, sample_no, is_forced_output);
                    }
                    else if (vcfRecord.is_snv())
                    {
                        posProcessor.insert_forced_output_pos(vcfRecord.pos - 1);
                    }
                }
                else if (INPUT_TYPE::KNOWN_VARIANTS == currentIndex)
                {
                    if (vcfRecord.is_indel())
                    {
                        processTrueIndelVariantRecord(opt.maxIndelSize, vcfRecord, posProcessor);
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
                if (INPUT_TYPE::EXCLUDE_REGION == currentIndex)
                {
                    const known_pos_range2 excludedRange(bedRecord.begin, bedRecord.end);
                    posProcessor.insertExcludedRegion(excludedRange);
                }
                else
                {
                    assert(false && "Unexpected hts index");
                }
            }
            else
            {
                assert(false && "Invalid input condition");
            }
        }
    }
    posProcessor.completeProcessing();
}
