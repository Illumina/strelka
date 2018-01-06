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

#include "starling_run.hh"
#include "starling_pos_processor.hh"
#include "starling_streams.hh"

#include "appstats/RunStatsManager.hh"
#include "blt_util/id_map.hh"
#include "blt_util/log.hh"
#include "common/Exceptions.hh"
#include "htsapi/bam_header_util.hh"
#include "htsapi/vcf_record_util.hh"
#include "starling_common/HtsMergeStreamerUtil.hh"
#include "starling_common/ploidy_util.hh"
#include "starling_common/starling_pos_processor_util.hh"



namespace INPUT_TYPE
{
enum index_t
{
    CANDIDATE_INDELS,
    FORCED_GT_VARIANTS,
    PLOIDY_REGION,
    NOCOMPRESS_REGION,
    CALL_REGION
};
}


/// \brief Parse sample names from ploidy VCF and match these to expected sample names
///
/// 1. Validate sample count/sample name match to input alignment files
/// 2. Construct a sample index translation map
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
    std::vector<unsigned>& sampleIndexToPloidyVcfSampleIndex)
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
            oss << "Repeated entry for sample name '" << vcfSampleName << "' in ploidy VCF file: '" << vcfStream.name() << "'";
            BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
        }

        vcfSampleSet.insert_key(vcfSampleName);
    }

    sampleIndexToPloidyVcfSampleIndex.clear();
    for (const auto& sampleName : sampleNames)
    {
        const auto maybeId(vcfSampleSet.get_optional_id(sampleName));
        if (not maybeId)
        {
            std::ostringstream oss;
            oss << "No entry for sample name '" << sampleName << "' in ploidy VCF file: '" << vcfStream.name() << "'";
            BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
        }
        else
        {
            sampleIndexToPloidyVcfSampleIndex.push_back(*maybeId);
        }
    }
}



static
void
callRegion(
    const starling_options& opt,
    const AnalysisRegionInfo& regionInfo,
    const starling_streams& fileStreams,
    const std::vector<unsigned>& sampleIndexToPloidyVcfSampleIndex,
    const unsigned ploidyVcfSampleCount,
    starling_read_counts& readCounts,
    reference_contig_segment& ref,
    HtsMergeStreamer& streamData,
    starling_pos_processor& posProcessor)
{
    using namespace illumina::common;

    const unsigned sampleCount(opt.alignFileOpt.alignmentFilenames.size());

    posProcessor.resetRegion(regionInfo.regionChrom, regionInfo.regionRange);
    streamData.resetRegion(regionInfo.streamerRegion.c_str());
    setRefSegment(opt, regionInfo.regionChrom, regionInfo.refRegionRange, ref);

    while (streamData.next())
    {
        const pos_t currentPos(streamData.getCurrentPos());
        const HTS_TYPE::index_t currentHtsType(streamData.getCurrentType());
        const unsigned currentIndex(streamData.getCurrentIndex());

        // wind posProcessor forward to position behind buffer head:
        posProcessor.set_head_pos(currentPos-1);

        if       (HTS_TYPE::BAM == currentHtsType)
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

            processInputReadAlignment(opt, ref, streamData.getCurrentBamStreamer(),
                                      streamData.getCurrentBam(), currentPos,
                                      readCounts, posProcessor, currentIndex);
        }
        else if (HTS_TYPE::VCF == currentHtsType)
        {
            assertExpectedVcfReference(ref, streamData.getCurrentVcfStreamer());
            const vcf_record& vcfRecord(streamData.getCurrentVcf());
            if     (INPUT_TYPE::CANDIDATE_INDELS == currentIndex)     // process candidate indels input from vcf file(s)
            {
                if (vcfRecord.is_indel())
                {
                    process_candidate_indel(opt.maxIndelSize, vcfRecord, posProcessor);
                }
                else
                {
                    log_os << "WARNING: candidate indel vcf variant record cannot be categorized as indel:\n";
                    streamData.getCurrentVcfStreamer().report_state(log_os);
                }
            }
            else if (INPUT_TYPE::FORCED_GT_VARIANTS == currentIndex)     // process forced genotype tests from vcf file(s)
            {
                if (vcfRecord.is_indel())
                {
                    static const unsigned sample_no(0);
                    static const bool is_forced_output(true);
                    process_candidate_indel(opt.maxIndelSize, vcfRecord, posProcessor, sample_no, is_forced_output);
                }
                else if (vcfRecord.is_snv() or vcfRecord.is_ref_site())
                {
                    posProcessor.insert_forced_output_pos(vcfRecord.pos - 1);
                }
                else
                {
                    std::ostringstream oss;
                    oss << "forcedGT vcf variant record cannot be categorized as SNV or indel:\n";
                    streamData.getCurrentVcfStreamer().report_state(oss);
                    BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
                }
            }
            else if (INPUT_TYPE::PLOIDY_REGION == currentIndex)
            {
                std::vector<unsigned> samplePloidy;
                known_pos_range2 ploidyRange;
                try
                {
                    parsePloidyFromVcf(ploidyVcfSampleCount, vcfRecord.line, ploidyRange, samplePloidy);
                }
                catch (...)
                {
                    log_os << "ERROR: Exception caught while parsing vcf ploidy record\n";
                    streamData.getCurrentVcfStreamer().report_state(log_os);
                    throw;
                }

                for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
                {
                    const unsigned ploidy(samplePloidy[sampleIndexToPloidyVcfSampleIndex[sampleIndex]]);
                    if ((ploidy == 0) || (ploidy == 1))
                    {
                        const bool retval(posProcessor.insert_ploidy_region(sampleIndex, ploidyRange, ploidy));
                        if (!retval)
                        {
                            std::ostringstream oss;
                            const auto& sampleName(fileStreams.getSampleNames()[sampleIndex]);
                            oss << "Ploidy vcf FORMAT/CN values conflict. Conflict detected in sample '"
                                << sampleName << "' at:\n";
                            streamData.getCurrentVcfStreamer().report_state(oss);
                            BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
                        }
                    }
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
            if (INPUT_TYPE::NOCOMPRESS_REGION == currentIndex)
            {
                known_pos_range2 range(bedRecord.begin,bedRecord.end);
                posProcessor.insert_nocompress_region(range);
            }
            else if (INPUT_TYPE::CALL_REGION == currentIndex)
            {
                known_pos_range2 range(bedRecord.begin,bedRecord.end);
                posProcessor.insertCallRegion(range);
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



void
starling_run(
    const prog_info& pinfo,
    const starling_options& opt)
{
    using namespace illumina::common;

    // ensure that this object is created first for runtime benchmark
    RunStatsManager statsManager(opt.segmentStatsFilename);

    opt.validate();

    const starling_deriv_options dopt(opt);
    starling_read_counts readCounts;
    reference_contig_segment ref;

    const unsigned sampleCount(opt.getSampleCount());

    ////////////////////////////////////////
    // setup streamData:
    //
    HtsMergeStreamer streamData(opt.referenceFilename);

    // additional data structures required in the region loop below, which are filled in as a side effect of
    // streamData initialization:
    std::vector<std::reference_wrapper<const bam_hdr_t>> bamHeaders;
    std::vector<std::string> sampleNames;
    unsigned ploidyVcfSampleCount(0);
    std::vector<unsigned> sampleIndexToPloidyVcfSampleIndex;

    {
        std::vector<unsigned> registrationIndices;
        for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
        {
            registrationIndices.push_back(sampleIndex);
        }
        bamHeaders = registerAlignments(opt.alignFileOpt.alignmentFilenames, registrationIndices, streamData);

        assert(not bamHeaders.empty());
        const bam_hdr_t& referenceHeader(bamHeaders.front());

        static const bool noRequireNormalized(false);
        registerVcfList(opt.input_candidate_indel_vcf, INPUT_TYPE::CANDIDATE_INDELS, referenceHeader, streamData,
                        noRequireNormalized);
        registerVcfList(opt.force_output_vcf, INPUT_TYPE::FORCED_GT_VARIANTS, referenceHeader, streamData);

        // initialize sampleNames from all bam headers (assuming 1 sample per bam for now)
        assert(bamHeaders.size() == sampleCount);
        for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
        {
            std::ostringstream defaultName;
            defaultName << "SAMPLE" << (sampleIndex+1);
            std::string sampleName(get_bam_header_sample_name(bamHeaders[sampleIndex], defaultName.str().c_str()));
            // remove spaces from sample name
            std::replace(sampleName.begin(), sampleName.end(), ' ', '_');
            sampleNames.push_back(sampleName);
        }

        if (!opt.ploidy_region_vcf.empty())
        {
            const vcf_streamer& vcfStream(streamData.registerVcf(opt.ploidy_region_vcf.c_str(), INPUT_TYPE::PLOIDY_REGION));
            vcfStream.validateBamHeaderChromSync(referenceHeader);

            mapVcfSampleIndices(vcfStream, sampleNames, sampleIndexToPloidyVcfSampleIndex);
            ploidyVcfSampleCount = vcfStream.getSampleCount();
        }

        if (!opt.gvcf.nocompress_region_bedfile.empty())
        {
            streamData.registerBed(opt.gvcf.nocompress_region_bedfile.c_str(), INPUT_TYPE::NOCOMPRESS_REGION);
        }

        if (! opt.callRegionsBedFilename.empty())
        {
            streamData.registerBed(opt.callRegionsBedFilename.c_str(), INPUT_TYPE::CALL_REGION);
        }
    }

    starling_streams fileStreams(opt, pinfo, bamHeaders, sampleNames);
    starling_pos_processor posProcessor(opt, dopt, ref, fileStreams, statsManager);

    const bam_hdr_t& referenceHeader(bamHeaders.front());
    const bam_header_info referenceHeaderInfo(referenceHeader);

    // parse and sanity check regions
    unsigned supplementalRegionBorderSize(opt.maxIndelSize);
    if (opt.isHaplotypingEnabled)
    {
        supplementalRegionBorderSize += ActiveRegionProcessor::MaxRefSpanToBypassAssembly;
    }
    const auto& referenceAlignmentFilename(opt.alignFileOpt.alignmentFilenames.front());
    std::vector<AnalysisRegionInfo> regionInfoList;
    getStrelkaAnalysisRegions(opt, referenceAlignmentFilename, referenceHeaderInfo, supplementalRegionBorderSize, regionInfoList);

    for (const auto& regionInfo : regionInfoList)
    {
        if (not opt.isUseCallRegions())
        {
            callRegion(opt, regionInfo, fileStreams, sampleIndexToPloidyVcfSampleIndex, ploidyVcfSampleCount,
                       readCounts, ref, streamData, posProcessor);
        }
        else
        {
            std::vector<known_pos_range2> subRegionRanges;
            getSubRegionsFromBedTrack(opt.callRegionsBedFilename, regionInfo.regionChrom, regionInfo.regionRange, subRegionRanges);

            for (const auto& subRegionRange : subRegionRanges)
            {
                AnalysisRegionInfo subRegionInfo;
                getStrelkaAnalysisRegionInfo(regionInfo.regionChrom, subRegionRange.begin_pos(), subRegionRange.end_pos(),
                                             supplementalRegionBorderSize, subRegionInfo);
                callRegion(opt, subRegionInfo, fileStreams, sampleIndexToPloidyVcfSampleIndex, ploidyVcfSampleCount,
                           readCounts, ref, streamData, posProcessor);
            }
        }
    }
    posProcessor.reset();
}
