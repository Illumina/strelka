// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#include "pedicure_pos_processor.hh"
#include "pedicure_streams.hh"
#include "pedicure_run.hh"

#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"
#include "htsapi/bam_header_util.hh"
#include "htsapi/bam_streamer.hh"
#include "htsapi/vcf_streamer.hh"
#include "starling_common/starling_input_stream_handler.hh"
#include "starling_common/starling_ref_seq.hh"
#include "starling_common/starling_pos_processor_util.hh"

#include <sstream>



void
pedicure_run(
    const prog_info& pinfo,
    const pedicure_options& opt)
{
    reference_contig_segment ref;
    get_starling_ref_seq(opt,ref);
    const pedicure_deriv_options dopt(opt,ref);

    const pos_range& rlimit(dopt.report_range_limit);

    assert(opt.bam_filename.empty());

    const std::string bam_region(get_starling_bam_region_string(opt,dopt));

    typedef std::shared_ptr<bam_streamer> stream_ptr;
    std::vector<stream_ptr> bamStreams;

    // setup all alignment data for main scan loop:
    for (const std::string& afile : opt.alignFileOpt.alignmentFilename)
    {
        stream_ptr tmp(new bam_streamer(afile.c_str(), bam_region.c_str()));
        bamStreams.push_back(tmp);
    }

    // check bam header compatibility:
    const unsigned bamCount(bamStreams.size());
    if (bamCount > 1)
    {
        /// TODO: provide a better error exception for failed bam header check:
        const bam_header_t* compareHeader(bamStreams[0]->get_header());
        for (unsigned bamIndex(1); bamIndex<bamCount; ++bamIndex)
        {
            const bam_header_t* indexHeader(bamStreams[bamIndex]->get_header());
            if (! check_header_compatibility(compareHeader,indexHeader))
            {
                log_os << "ERROR: incompatible bam headers between files:\n"
                       << "\t" << opt.alignFileOpt.alignmentFilename[0] << "\n"
                       << "\t" << opt.alignFileOpt.alignmentFilename[bamIndex] << "\n";
                exit(EXIT_FAILURE);
            }
        }
    }

    const int32_t tid(bamStreams[0]->target_name_to_id(opt.bam_seq_name.c_str()));
    if (tid < 0)
    {
        std::ostringstream oss;
        oss << "ERROR: seq_name: '" << opt.bam_seq_name << "' is not found in the header of BAM file: '" <<  opt.alignFileOpt.alignmentFilename[0] << "'\n";
        throw blt_exception(oss.str().c_str());
    }

    // We make the assumption that all alignment files have the
    // same set of reference chromosomes (and thus tid matches for the
    // same chrom label in the binary records). Check this constraint
    // here:
    for (unsigned bamIndex(1); bamIndex<bamCount; ++bamIndex)
    {
        const int32_t other_tid(bamStreams[bamIndex]->target_name_to_id(opt.bam_seq_name.c_str()));
        if (tid != other_tid)
        {
            throw blt_exception("ERROR: sample BAM files have mis-matched reference sequence dictionaries.\n");
        }
    }

    const PedicureSampleSetSummary ssi(opt);
    const bam_header_t* const header(bamStreams[0]->get_header());
    pedicure_streams streams(opt, dopt, pinfo, header, ssi);
    pedicure_pos_processor sppr(opt,dopt,ref,streams);
    starling_read_counts brc;

    starling_input_stream_data sdata;
    for (unsigned bamIndex(0); bamIndex<bamCount; ++bamIndex)
    {
        sdata.register_reads(*bamStreams[bamIndex],bamIndex);
    }

    // hold zero-to-many vcf streams open:
    typedef std::shared_ptr<vcf_streamer> vcf_ptr;
    std::vector<vcf_ptr> indel_stream;

    for (const auto& vcf_filename : opt.input_candidate_indel_vcf)
    {
        indel_stream.push_back(vcf_ptr(new vcf_streamer(vcf_filename.c_str(),
                                                        bam_region.c_str(),header)));
        sdata.register_indels(*(indel_stream.back()));
    }

    std::vector<vcf_ptr> foutput_stream;

    for (const auto& vcf_filename : opt.force_output_vcf)
    {
        foutput_stream.push_back(vcf_ptr(new vcf_streamer(vcf_filename.c_str(),
                                                          bam_region.c_str(),header)));
        sdata.register_forced_output(*(foutput_stream.back()));
    }

    starling_input_stream_handler sinput(sdata);

    while (sinput.next())
    {
        const input_record_info current(sinput.get_current());

        // If we're past the end of rlimit range then we're done.
        //   Note that some additional padding is allowed for off
        //   range indels which might influence results within the
        //   report range:
        //
        if (rlimit.is_end_pos && (current.pos >= (rlimit.end_pos+static_cast<pos_t>(opt.max_indel_size)))) break;

        // wind sppr forward to position behind buffer head:
        sppr.set_head_pos(sinput.get_head_pos()-1);

        if       (current.itype == INPUT_TYPE::READ)   // handle reads from the primary mapper (as opposed to the local assembler)
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

            const bam_streamer& readStream(*bamStreams[current.sample_no]);
            const bam_record& read(*(readStream.get_record_ptr()));

            process_genomic_read(opt,ref,readStream,read,current.pos,
                                 rlimit.begin_pos,brc,sppr,current.sample_no);
        }
        else if (current.itype == INPUT_TYPE::INDEL)     // process candidate indels input from vcf file(s)
        {
            const vcf_record& vcf_indel(*(indel_stream[current.get_order()]->get_record_ptr()));
            process_candidate_indel(opt.max_indel_size, vcf_indel, sppr);

        }
        else if (current.itype == INPUT_TYPE::FORCED_OUTPUT)     // process forced genotype tests from vcf file(s)
        {
            const vcf_record& vcf_variant(*(foutput_stream[current.get_order()]->get_record_ptr()));
            if (vcf_variant.is_indel())
            {
                static const unsigned sample_no(0);
                static const bool is_forced_output(true);
                process_candidate_indel(opt.max_indel_size, vcf_variant,sppr,sample_no,is_forced_output);
            }
            else if (vcf_variant.is_snv())
            {
                sppr.insert_forced_output_pos(vcf_variant.pos-1);
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
