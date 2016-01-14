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

#include "snoise_run.hh"
#include "snoise_pos_processor.hh"
#include "snoise_streams.hh"

#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"
#include "htsapi/bam_streamer.hh"
#include "htsapi/vcf_streamer.hh"
#include "starling_common/starling_base_shared.hh"
#include "starling_common/starling_input_stream_handler.hh"
#include "starling_common/starling_ref_seq.hh"
#include "starling_common/starling_pos_processor_util.hh"
#include <sstream>



void
snoise_run(
    const prog_info& pinfo,
    const snoise_options& opt)
{
    reference_contig_segment ref;
    get_starling_ref_seq(opt,ref);

    const starling_base_deriv_options dopt(opt,ref);
    const pos_range& rlimit(dopt.report_range_limit);

    assert(! opt.bam_filename.empty());

    const std::string bam_region(get_starling_bam_region_string(opt,dopt));
    bam_streamer read_stream(opt.bam_filename.c_str(),bam_region.c_str());

    const int32_t tid(read_stream.target_name_to_id(opt.bam_seq_name.c_str()));
    if (tid < 0)
    {
        std::ostringstream oss;
        oss << "ERROR: seq_name: '" << opt.bam_seq_name << "' is not found in the header of BAM file: '" << opt.bam_filename << "'\n";
        throw blt_exception(oss.str().c_str());
    }

    SampleSetSummary ssi;
    snoise_streams streams(opt,pinfo,read_stream.get_header(),ssi);

    snoise_pos_processor sppr(opt,dopt,ref,streams);
    starling_read_counts brc;

    starling_input_stream_data sdata;
    sdata.register_reads(read_stream);

    // hold zero-to-many vcf streams open:
    typedef std::shared_ptr<vcf_streamer> vcf_ptr;
    std::vector<vcf_ptr> indel_stream;

    for (const auto& vcf_filename : opt.input_candidate_indel_vcf)
    {
        indel_stream.push_back(vcf_ptr(new vcf_streamer(vcf_filename.c_str(),
                                                        bam_region.c_str(),read_stream.get_header())));
        sdata.register_indels(*(indel_stream.back()));
    }

    std::vector<vcf_ptr> foutput_stream;

    for (const auto& vcf_filename : opt.force_output_vcf)
    {
        foutput_stream.push_back(vcf_ptr(new vcf_streamer(vcf_filename.c_str(),
                                                          bam_region.c_str(),read_stream.get_header())));
        sdata.register_forced_output(*(foutput_stream.back()));
    }

    starling_input_stream_handler sinput(sdata);

    while (sinput.next())
    {
        const input_record_info current(sinput.get_current());

        // Process finishes at the the end of rlimit range. Note that
        // some additional padding is allowed for off-range indels
        // which might influence results within rlimit:
        //
        if (rlimit.is_end_pos && (current.pos >= (rlimit.end_pos+static_cast<pos_t>(opt.max_indel_size)))) break;

        // wind sppr forward to position behind buffer head:
        sppr.set_head_pos(sinput.get_head_pos()-1);

        if       (current.itype == INPUT_TYPE::READ)
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

            const bam_record& read(*(read_stream.get_record_ptr()));

            process_genomic_read(opt,ref,read_stream,read,current.pos,rlimit.begin_pos,brc,sppr);

        }
        else if (current.itype == INPUT_TYPE::INDEL)     // process candidate indels input from vcf file(s)
        {
            const vcf_record& vcf_indel(*(indel_stream[current.get_order()]->get_record_ptr()));
            process_candidate_indel(opt.max_indel_size, vcf_indel,sppr);

        }
        else if (current.itype == INPUT_TYPE::FORCED_OUTPUT)     // process forced genotype tests from vcf file(s)
        {
            const vcf_record& vcf_variant(*(foutput_stream[current.get_order()]->get_record_ptr()));
            if       (vcf_variant.is_indel())
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

    //    brc.report(client_io.report_os());
}

