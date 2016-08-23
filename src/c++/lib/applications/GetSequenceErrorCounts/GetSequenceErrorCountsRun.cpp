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

#include "GetSequenceErrorCountsRun.hh"
#include "SequenceErrorCountsPosProcessor.hh"

#include "appstats/RunStatsManager.hh"
#include "blt_util/log.hh"
#include "common/Exceptions.hh"
#include "htsapi/align_path_bam_util.hh"
#include "starling_common/HtsMergeStreamerUtil.hh"
#include "starling_common/ploidy_util.hh"
#include "starling_common/starling_ref_seq.hh"
#include "starling_common/starling_pos_processor_util.hh"

#include <sstream>
#include <algorithm>


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


/// adds truth variant input from a vcf to indel synchronizer
///
static
void
processTrueIndelVariantRecord(
    const unsigned max_indel_size,
    const vcf_record& vcfr,
    SequenceErrorCountsPosProcessor& sppr)
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

            sppr.insert_indel(obs,sample_no);
            sppr.addKnownVariant(vcfr);
        }
    }
}



void
getSequenceErrorCountsRun(
    const prog_info& pinfo,
    const SequenceErrorCountsOptions& opt)
{
    opt.validate();

    RunStatsManager segmentStatMan(opt.segmentStatsFilename);

    reference_contig_segment ref;
    get_starling_ref_seq(opt,ref);

    const SequenceErrorCountsDerivOptions dopt(opt,ref);
    const pos_range& rlimit(dopt.report_range_limit);

    assert(! opt.bam_filename.empty());

    const std::string bam_region(get_starling_bam_region_string(opt,dopt));
    HtsMergeStreamer streamData(bam_region.c_str());
    const bam_streamer& readStream(streamData.registerBam(opt.bam_filename.c_str()));
    const bam_hdr_t& readHeader(readStream.get_header());

    const int32_t tid(readStream.target_name_to_id(opt.bam_seq_name.c_str()));
    if (tid < 0)
    {
        using namespace illumina::common;

        std::ostringstream oss;
        oss << "ERROR: seq_name: '" << opt.bam_seq_name << "' is not found in the header of BAM/CRAM file: '" << opt.bam_filename << "'\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    SampleSetSummary ssi;
    SequenceErrorCountsStreams client_io(opt,pinfo,readHeader,ssi);

    SequenceErrorCountsPosProcessor sppr(opt,dopt,ref,client_io);
    starling_read_counts brc;

    registerVcfList(opt.input_candidate_indel_vcf, INPUT_TYPE::CANDIDATE_INDELS, readHeader, streamData);
    registerVcfList(opt.force_output_vcf, INPUT_TYPE::FORCED_GT_VARIANTS, readHeader, streamData);

    if (! opt.knownVariantsFile.empty())
    {
        const vcf_streamer& vcfStream(streamData.registerVcf(opt.knownVariantsFile.c_str(), INPUT_TYPE::KNOWN_VARIANTS));
        vcfStream.validateBamHeaderChromSync(readHeader);
    }

    for (const std::string& excludeRegionFilename : opt.excludedRegionsFileList)
    {
        streamData.registerBed(excludeRegionFilename.c_str(), INPUT_TYPE::EXCLUDE_REGION);
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

            const bam_record& read(streamData.getCurrentBam());

            // special test used in error counting only -- this isn't the ideal place for this:
            {
                using namespace ALIGNPATH;
                path_t apath;
                bam_cigar_to_apath(read.raw_cigar(),read.n_cigar(),apath);

                if (apath_indel_count(apath) > 2) continue;
            }

            processInputReadAlignment(opt, ref, streamData.getCurrentBamStreamer(),
                                      read, currentPos, rlimit.begin_pos, brc, sppr);
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
            else if (INPUT_TYPE::KNOWN_VARIANTS == currentIndex)
            {
                if (vcfRecord.is_indel())
                {
                    processTrueIndelVariantRecord(opt.max_indel_size, vcfRecord, sppr);
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
            if     (INPUT_TYPE::EXCLUDE_REGION == currentIndex)
            {
                const known_pos_range2 excludedRange(bedRecord.begin,bedRecord.end);
                sppr.insertExcludedRegion(excludedRange);
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
