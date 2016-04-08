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

#include "strelka_pos_processor.hh"
#include "strelka_streams.hh"
#include "strelka_run.hh"

#include "appstats/RunStatsManager.hh"
#include "blt_util/log.hh"
#include "htsapi/bam_header_util.hh"
#include "starling_common/HtsMergeStreamerUtil.hh"
#include "starling_common/starling_ref_seq.hh"
#include "starling_common/starling_pos_processor_util.hh"

#include <sstream>



namespace INPUT_TYPE
{
enum index_t
{
    CANDIDATE_INDELS,
    FORCED_GT_VARIANTS,
    NOISE_VARIANTS,
    PLOIDY_REGION,
    NOCOMPRESS_REGION
};
}



void
strelka_run(
    const prog_info& pinfo,
    const strelka_options& opt)
{
    opt.validate();

    RunStatsManager segmentStatMan(opt.segmentStatsFilename);

    reference_contig_segment ref;
    get_starling_ref_seq(opt,ref);
    const strelka_deriv_options dopt(opt,ref);

    const pos_range& rlimit(dopt.report_range_limit);

    assert(! opt.tumor_bam_filename.empty());

    const std::string bam_region(get_starling_bam_region_string(opt,dopt));
    HtsMergeStreamer streamData(bam_region.c_str());
    const bam_streamer& tumorReadStream(streamData.registerBam(opt.tumor_bam_filename.c_str(),STRELKA_SAMPLE_TYPE::TUMOR));
    const bam_hdr_t& tumorReadHeader(tumorReadStream.get_header());

    std::unique_ptr<bam_streamer> normal_read_stream_ptr;

    if (! opt.bam_filename.empty())
    {
        const bam_streamer& normalReadStream(streamData.registerBam(opt.bam_filename.c_str(),STRELKA_SAMPLE_TYPE::NORMAL));

        if (! check_header_compatibility(normalReadStream.get_header(),tumorReadHeader))
        {
            std::ostringstream oss;
            oss << "ERROR: Normal and tumor BAM/CRAM files have incompatible headers.\n";
            oss << "\tnormal_bam_file:\t'" << opt.bam_filename << "'\n";
            oss << "\ttumor_bam_file:\t'" << opt.tumor_bam_filename << "'\n";
            throw blt_exception(oss.str().c_str());
        }
    }

    const int32_t tid(tumorReadStream.target_name_to_id(opt.bam_seq_name.c_str()));
    if (tid < 0)
    {
        std::ostringstream oss;
        oss << "ERROR: seq_name: '" << opt.bam_seq_name << "' is not found in the header of BAM/CRAM file: '" << opt.tumor_bam_filename << "'\n";
        throw blt_exception(oss.str().c_str());
    }

    const StrelkaSampleSetSummary ssi;
    strelka_streams client_io(opt, dopt, pinfo, tumorReadHeader,ssi);
    strelka_pos_processor sppr(opt,dopt,ref,client_io);
    starling_read_counts brc;

    registerVcfList(opt.input_candidate_indel_vcf, INPUT_TYPE::CANDIDATE_INDELS, tumorReadHeader, streamData);
    registerVcfList(opt.force_output_vcf, INPUT_TYPE::FORCED_GT_VARIANTS, tumorReadHeader, streamData);
    registerVcfList(opt.noise_vcf, INPUT_TYPE::NOISE_VARIANTS, tumorReadHeader, streamData);

    while (streamData.next())
    {
        const pos_t currentPos(streamData.getCurrentPos());
        const HTS_TYPE::index_t currentHtsType(streamData.getCurrentType());
        const unsigned currentIndex(streamData.getCurrentIndex());

        // If we're past the end of rlimit range then we're done.
        //   Note that some additional padding is allowed for off
        //   range indels which might influence results within the
        //   report range:
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
            //if((current_pos+MAX_READ_SIZE+MAX_INDEL_SIZE) <= rlimit.begin_pos) continue;
            process_genomic_read(opt,ref,streamData.getCurrentBamStreamer(),
                                 streamData.getCurrentBam(),currentPos,
                                 rlimit.begin_pos,brc,sppr,currentIndex);
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
            else if (INPUT_TYPE::NOISE_VARIANTS == currentIndex)
            {
                if (vcfRecord.is_snv())
                {
                    SiteNoise sn;
                    set_noise_from_vcf(vcfRecord.line,sn);
                    sppr.insert_noise_pos(vcfRecord.pos-1,sn);
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

    //    brc.report(client_io.report_os());
}
