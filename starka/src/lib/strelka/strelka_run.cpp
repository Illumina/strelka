// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

/// \file

/// \author Chris Saunders
///

#include "strelka_pos_processor.hh"
#include "strelka_streams.hh"

#include "blt_util/bam_streamer.hh"
#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"
#include "blt_util/vcf_streamer.hh"
#include "starling_common/grouper_contig_util.hh"
#include "starling_common/starling_input_stream_handler.hh"
#include "starling_common/starling_ref_seq.hh"
#include "starling_common/starling_pos_processor_contig_util.hh"
#include "starling_common/starling_pos_processor_util.hh"
#include "strelka/strelka_info.hh"
#include "strelka/strelka_run.hh"

#include "boost/shared_ptr.hpp"

#include <sstream>


namespace {
const prog_info& pinfo(strelka_info::get());
}



void
strelka_run(const strelka_options& opt) {

    reference_contig_segment ref;
    get_starling_ref_seq(opt,ref);

    const strelka_deriv_options dopt(opt,ref);

    const pos_range& rlimit(dopt.report_range_limit);

    contig_data_manager normal_cdm(opt.indel_contig_filename,
                                   opt.indel_contig_read_filename);

    contig_data_manager tumor_cdm(opt.tumor_indel_contig_filename,
                                  opt.tumor_indel_contig_read_filename);

    assert(! opt.bam_filename.empty());
    assert(! opt.tumor_bam_filename.empty());

    const std::string bam_region(get_starling_bam_region_string(opt,dopt));
    bam_streamer normal_read_stream(opt.bam_filename.c_str(),bam_region.c_str());
    bam_streamer tumor_read_stream(opt.tumor_bam_filename.c_str(),bam_region.c_str());

    // check for header consistency:
    if (! check_header_compatibility(normal_read_stream.get_header(),tumor_read_stream.get_header())) {
        std::ostringstream oss;
        oss << "ERROR: Normal and tumor BAM files have incompatible headers.\n";
        oss << "\tnormal_bam_file:\t'" << opt.bam_filename << "'\n";
        oss << "\ttumor_bam_file:\t'" << opt.tumor_bam_filename << "'\n";
        throw blt_exception(oss.str().c_str());
    }

    const int32_t tid(normal_read_stream.target_name_to_id(opt.bam_seq_name.c_str()));
    if (tid < 0) {
        std::ostringstream oss;
        oss << "ERROR: seq_name: '" << opt.bam_seq_name << "' is not found in the header of BAM file: '" << opt.bam_filename << "'\n";
        throw blt_exception(oss.str().c_str());
    }

    // We make the assumption that the normal and tumor files have the
    // same set of reference chromosomes (and thus tid matches for the
    // same chrom label in the binary records). Check this constraint
    // here:
    {
        const int32_t tumor_tid(tumor_read_stream.target_name_to_id(opt.bam_seq_name.c_str()));
        if (tid != tumor_tid) {
            throw blt_exception("ERROR: tumor and normal BAM files have mis-matched reference sequence dictionaries.\n");
        }
    }

    // Provide a temporary bam record for contig reads to write key
    // information into, and initialize this record with values that
    // will be fixed for this run:
    bam_record tmp_key_br;
    tmp_key_br.set_target_id(tid);

    strelka_streams client_io(opt,pinfo,normal_read_stream.get_header());
    strelka_pos_processor sppr(opt,dopt,ref,client_io);
    starling_read_counts brc;

    starling_input_stream_data sdata;
    sdata.register_reads(normal_read_stream,STRELKA_SAMPLE_TYPE::NORMAL);
    sdata.register_reads(tumor_read_stream,STRELKA_SAMPLE_TYPE::TUMOR);
    sdata.register_contigs(normal_cdm.creader(),STRELKA_SAMPLE_TYPE::NORMAL);
    sdata.register_contigs(tumor_cdm.creader(),STRELKA_SAMPLE_TYPE::TUMOR);

    // hold zero-to-many vcf streams open in indel_streams:
    typedef boost::shared_ptr<vcf_streamer> vcf_ptr;
    std::vector<vcf_ptr> indel_stream;
    for (unsigned i(0); i<opt.input_candidate_indel_vcf.size(); ++i) {
        indel_stream.push_back(vcf_ptr(new vcf_streamer(opt.input_candidate_indel_vcf[i].c_str(),
                                                        bam_region.c_str(),normal_read_stream.get_header())));
        sdata.register_indels(*(indel_stream.back()));
    }

    starling_input_stream_handler sinput(sdata);

    while (sinput.next()) {
        const input_record_info current(sinput.get_current());

        // If we're past the end of rlimit range then we're done.
        //   Note that some additional padding is allowed for off
        //   range indels which might influence results within the
        //   report range:
        //
        if (rlimit.is_end_pos && (current.pos >= (rlimit.end_pos+static_cast<pos_t>(opt.max_indel_size)))) break;

        // wind sppr forward to position behind buffer head:
        sppr.set_head_pos(sinput.get_head_pos()-1);

        if       (current.itype == INPUT_TYPE::READ) { // handle reads from the primary mapper (as opposed to the local assembler)

            // Remove the filter below because it's not valid for
            // RNA-Seq case, reads should be selected for the report
            // range by the bam reading functions
            //
            // /// get potential bounds of the read based only on current_pos:
            // const known_pos_range any_read_bounds(current_pos-max_indel_size,current_pos+MAX_READ_SIZE+max_indel_size);
            // if( sppr.is_range_outside_report_influence_zone(any_read_bounds) ) continue;

            // Approximate begin range filter: (removed for RNA-Seq)
            //if((current_pos+MAX_READ_SIZE+MAX_INDEL_SIZE) <= rlimit.begin_pos) continue;
            const bam_streamer* streamptr(NULL);
            if        (current.sample_no == STRELKA_SAMPLE_TYPE::NORMAL) {
                streamptr = &normal_read_stream;
            } else if (current.sample_no == STRELKA_SAMPLE_TYPE::TUMOR) {
                streamptr = &tumor_read_stream;
            } else {
                log_os << "ERROR: unrecognized sample_no: " << current.sample_no << "\n";
                exit(EXIT_FAILURE);
            }
            const bam_streamer& read_stream(*streamptr);
            const bam_record& read(*(read_stream.get_record_ptr()));

            process_genomic_read(opt,ref,read_stream,read,current.pos,
                                 rlimit.begin_pos,brc,sppr,current.sample_no);

        } else if (current.itype == INPUT_TYPE::CONTIG) { // process local-assembly contig and its reads

            contig_data_manager* cdmp(NULL);
            if        (current.sample_no == STRELKA_SAMPLE_TYPE::NORMAL) {
                cdmp = &normal_cdm;
            } else if (current.sample_no == STRELKA_SAMPLE_TYPE::TUMOR) {
                cdmp = &tumor_cdm;
            } else {
                log_os << "ERROR: unrecognized sample_no: " << current.sample_no << "\n";
                exit(EXIT_FAILURE);
            }

            const grouper_contig& ctg(cdmp->creader().get_contig());

            const char* sample_label(STRELKA_SAMPLE_TYPE::get_label(current.sample_no));

            if (! test_contig_usability(opt,ctg,sppr,sample_label)) continue;

            process_contig(opt,ref,ctg,sppr,current.sample_no,sample_label);

            process_contig_reads(ctg,opt.max_indel_size,cdmp->contig_read_exr(),sppr,tmp_key_br,current.sample_no);

        } else if (current.itype == INPUT_TYPE::INDEL) { // process candidate indels input from vcf file(s)
            const vcf_record& vcf_indel(*(indel_stream[current.get_order()]->get_record_ptr()));
            process_candidate_indel(vcf_indel,sppr);

        } else {
            log_os << "ERROR: invalid input condition.\n";
            exit(EXIT_FAILURE);
        }
    }

    sppr.reset();

    //    brc.report(client_io.report_os());
}

