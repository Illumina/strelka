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
///
/// \author Chris Saunders
///
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the begining of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///


#include "candidate_alignment.hh"
#include "starling_pos_processor_indel_util.hh"
#include "starling_read_util.hh"

#include "blt_util/log.hh"
#include "starling_common/grouper_contig_util.hh"
#include "starling_common/starling_pos_processor_contig_util.hh"

#include "boost/foreach.hpp"
#include "boost/scoped_array.hpp"



static
void
print_contig_alignment(const grouper_contig& ctg,
                       std::ostream& os) {
    std::string cigar;
    apath_to_cigar(ctg.path,cigar);
    os << "\tContig alignment: '" << cigar << "'\n";
}



bool
test_contig_usability(const starling_options& opt,
                      const grouper_contig& ctg,
                      const starling_pos_processor_base& sppr,
                      const char* sample_label) {

    // skip contigs which are too large:
    if(! ctg.is_usable) return false;

    {   // skip contigs which will not affect results in the report range:
        const pos_t contig_span(apath_ref_length(ctg.path));
        const pos_t contig_end_pos(ctg.pos+contig_span);
        const known_pos_range contig_bounds(ctg.pos,contig_end_pos);
        const bool is_skip_contig(sppr.is_range_outside_report_influence_zone(contig_bounds));
        if(is_skip_contig) return false;
    }

    // skip contigs which contain adjacent insertion/deletion events:
    static const bool is_filter_swap(false);
    if(is_filter_swap && is_seq_swap(ctg.path)) {
        log_os << "WARNING: Indel contig: '" << ctg.id << "'";
        if(NULL != sample_label) log_os << " from sample '" << sample_label << "'";
        log_os << " contains adjacent insertion/deletion event. Skipping...\n";
        return false;
    }

    // skip contigs without a continuous segment of non-indel bases of
    // at least length N at each edge of the alignment:
    if(opt.min_contig_edge_alignment) {
        bool is_skip(false);
        const unsigned aps(ctg.path.size());
        const std::pair<unsigned,unsigned> res(get_match_edge_segments(ctg.path));
        assert((aps != res.first) && (aps != res.second));
        if((ctg.path[res.first].length < opt.min_contig_edge_alignment) ||
           (ctg.path[res.second].length < opt.min_contig_edge_alignment)) {
            is_skip=true;
        }
        if(is_skip) {
            if(opt.verbosity >= LOG_LEVEL::ALLWARN) {
                log_os << "WARNING: Indel contig: '" << ctg.id << "'";
                if(NULL != sample_label) log_os << " from sample '" << sample_label << "'";
                log_os << " contains an edge match segment shorter than "
                       << opt.min_contig_edge_alignment << " bases. Skipping...\n";
                print_contig_alignment(ctg,log_os);
            }
            return false;
        }
    }

    // skip contigs without a match segment of length at least N:
    if(opt.min_contig_contiguous_match) {
        using namespace ALIGNPATH;

        bool is_skip(true);
        BOOST_FOREACH(const path_segment& ps, ctg.path) {
            if((MATCH == ps.type) &&
               (opt.min_contig_contiguous_match <= ps.length)) {
                is_skip=false;
                break;
            }
        }
        if(is_skip) {
            if(opt.verbosity >= LOG_LEVEL::ALLWARN) {
                log_os << "WARNING: Indel contig: '" << ctg.id << "'";
                if(NULL != sample_label) log_os << " from sample '" << sample_label << "'";
                log_os << " does not contain a contiguous match segment of at least "
                       << opt.min_contig_contiguous_match << " bases. Skipping...\n";
                print_contig_alignment(ctg,log_os);
            }
            return false;
        }
    }

    // skip contigs with very short open-ended insertions at the end,
    // GROUPER currently generates many of these, so put this into the
    // extended warning set:
    if(opt.min_contig_open_end_support) {
        bool is_skip(false);
        const unsigned clead(apath_insert_lead_size(ctg.path));
        const unsigned ctrail(apath_insert_trail_size(ctg.path));
        if(clead && (clead<opt.min_contig_open_end_support)) is_skip=true;
        if(ctrail && (ctrail<opt.min_contig_open_end_support)) is_skip=true;
        if(is_skip) {
            if(opt.verbosity >= LOG_LEVEL::ALLWARN) {
                log_os << "WARNING: Indel contig: '" << ctg.id << "'";
                if(NULL != sample_label) log_os << " from sample '" << sample_label << "'";
                log_os << " contains short open-end. Skipping...\n";
                print_contig_alignment(ctg,log_os);
            }
            return false;
        }
    }

    return true;
}



/// insert an alignment record in export format into the
/// starling_processor
///
/// returns true if the read is accepted into the buffer
/// (read can fail a number of quality checks -- such as
/// being located too far away from other alignments of the
/// same read or having an indel that is too large
///
static
std::pair<bool,align_id_t>
add_export_read_to_buffer(starling_pos_processor_base& sppr,
                          const export_line_parser& exl,
                          const char* fwd_strand_read,
                          const unsigned read_size,
                          const alignment& al,
                          const READ_ALIGN::index_t rat,
                          const char* chrom_name,
                          const MAPLEVEL::index_t maplev,
                          bam_record& tmp_key_br,
                          const unsigned sample_no,
                          const align_id_t contig_no,
                          const indel_set_t* contig_indels_ptr) {

    // clear certain tmp_key_br fields, note that this temporary is
    // only used for export conversion so we shouldn't have to worry
    // about aux data, but we take care of it here just be be
    // thorough.
    //
    {
        bam1_t* brp(tmp_key_br.get_data());
        bam1_t& br(*brp);
        bam1_core_t& brc(br.core);
        brc.flag = 0;
        brc.pos = -1;
        brc.mpos = -1;
        // always set to the same chrom: brc.tid = -1;
        brc.mtid = -1;
        brc.isize = 0;

        br.data_len -= br.l_aux;
        br.l_aux = 0;
    }

    // set qname:
    std::string key_name;
    get_read_key_from_export_line(exl,key_name);
    tmp_key_br.set_qname(key_name.c_str());

    // set paired bit:
    if(exl.is_partner_strand() != tmp_key_br.is_paired()) {
        tmp_key_br.toggle_is_paired();
    }

    // set readno:
    const int readnum(exl.read_number());
    assert((readnum == 1) || (readnum == 2));
    if(readnum==1) {
        if(! tmp_key_br.is_first()) tmp_key_br.toggle_is_first();
        if(tmp_key_br.is_second()) tmp_key_br.toggle_is_second();
    } else {
        if(tmp_key_br.is_first()) tmp_key_br.toggle_is_first();
        if(! tmp_key_br.is_second()) tmp_key_br.toggle_is_second();
    }

    // set strand:
    if(al.is_fwd_strand!=tmp_key_br.is_fwd_strand()) {
        tmp_key_br.toggle_is_fwd_strand();
    }

    // set read and quality:
    const char* const raw_quality(exl.quality());
    boost::scoped_array<uint8_t> fwd_quality(new uint8_t[read_size]);
    if(al.is_fwd_strand) {
        for(unsigned i(0); i<read_size; ++i) {
            fwd_quality[i] = raw_quality[i]-64;
        }
    } else {
        for(unsigned i(0); i<read_size; ++i) {
            fwd_quality[i] = raw_quality[read_size-(i+1)]-64;
        }
    }

    tmp_key_br.set_readqual(fwd_strand_read,fwd_quality.get());

    // set mate strand:
    bool is_mate_fwd(true);
    if(exl.is_partner_strand()) {
        const char dir(exl.partner_strand());
        if(dir=='N') {
            is_mate_fwd=(exl.match_strand()!='F');
        } else {
            is_mate_fwd=(dir=='F');
        }
    }
    if(is_mate_fwd != tmp_key_br.is_mate_fwd_strand()) {
        tmp_key_br.toggle_is_mate_fwd_strand();
    }

    // \TODO add more settings for MPOS and ISIZE

    return sppr.insert_read(tmp_key_br,
                            al,
                            rat,
                            chrom_name,
                            maplev,
                            sample_no,
                            contig_no,
                            contig_indels_ptr);
}



// advance to first read aligning to GROUPER contig_id:
//
static
void
find_first_read_for_contig(const std::string& contig_id,
                           export_stream_reader& exr) {

    do {
        if((exr.exline()) &&
           (contig_id == exr.exline()->match_chromosome())) return;
    } while(exr.next());

    log_os << "ERROR: expected at least one read for GROUPER contig: " << contig_id << "\n";
    log_os << "       in GROUPER read export stream: " <<  exr.name() << "\n";
    exit(EXIT_FAILURE);
}



/// Mark any edge insertions in alignment as belonging to breakpoints.
///
/// Note this is only intended for assembled contig alignments.
///
static
void
set_candidate_edge_insert_bp(candidate_alignment& cal) {

    using namespace ALIGNPATH;

    const path_t& path(cal.al.path);

    pos_t ref_head_pos(cal.al.pos);
    const unsigned as(path.size());
    for(unsigned i(0); i<as; ++i) {
        const path_segment& ps(path[i]);
        const bool is_edge_segment((i==0) || ((i+1)==as));
        if((ps.type==INSERT) || (ps.type==DELETE)) {
            if(is_edge_segment) {
                assert(ps.type==INSERT);
                if(i==0) {
                    assert(cal.leading_indel_key.type == INDEL::NONE);
                    cal.leading_indel_key.type = INDEL::BP_RIGHT;
                    cal.leading_indel_key.pos = ref_head_pos;
                    cal.leading_indel_key.length = ps.length;
                } else {
                    assert(cal.trailing_indel_key.type == INDEL::NONE);
                    cal.trailing_indel_key.type = INDEL::BP_LEFT;
                    cal.trailing_indel_key.pos = ref_head_pos;
                    cal.trailing_indel_key.length = ps.length;
                }
            }
        } else {
            // routine assumes no clipping:
            assert(ps.type==MATCH);
        }
        if(is_segment_type_ref_length(ps.type)) ref_head_pos += ps.length;
    }
}


namespace {

// simple object to share code between the genomic and
// contig read parse routines
//
struct export_read_parse  : private boost::noncopyable {

    export_read_parse(const export_stream_reader& exr)
        : _fwd_strand_read_store(NULL)
    {
        const export_line_parser& exl(*(exr.exline()));

        read=exl.read();
        assert(read);

        read_size=strlen(read);
        assert(read_size>0);

        if(read_size > STARLING_MAX_READ_SIZE) {
            log_os << "ERROR: maximum read size (" << STARLING_MAX_READ_SIZE << ") exceeded in export line.\n";
            exr.report_state(log_os);
            exit(EXIT_FAILURE);
        }

        if(! is_valid_seq(read)) {
            log_os << "ERROR: unsupported base(s) in read sequence.\n";
            exr.report_state(log_os);
            exit(EXIT_FAILURE);
        }

        is_fwd_strand=(exl.match_strand() == 'F');

        fwd_strand_read=read;
        if(! is_fwd_strand) {
            _fwd_strand_read_store=new char[read_size+1];
            reverseCompCopy(read,read+read_size,_fwd_strand_read_store);
            fwd_strand_read=_fwd_strand_read_store;
        }
    }

    ~export_read_parse() {
        if(NULL!=_fwd_strand_read_store) delete [] _fwd_strand_read_store;
    }

    const char* read;
    unsigned read_size;
    const char* fwd_strand_read;
    bool is_fwd_strand;

private:
    char* _fwd_strand_read_store;
};



} // namespace


// TODO -- do I have to say it?
align_id_t contig_no(0);



void
process_contig_reads(const grouper_contig& ctg,
                     const unsigned max_indel_size,
                     export_stream_reader& exr,
                     starling_pos_processor_base& sppr,
                     bam_record& tmp_key_br,
                     const unsigned sample_no) {

    // advance to first read aligning to contig id:
    find_first_read_for_contig(ctg.id,exr);

    // get contig indel set:
    indel_set_t contig_indels;
    {
        candidate_alignment ctg_cal;
        ctg_cal.al=static_cast<alignment>(ctg);
        set_candidate_edge_insert_bp(ctg_cal);
        get_alignment_indels(ctg_cal,max_indel_size,contig_indels);
    }

    do {
        const export_line_parser& exl(*(exr.exline()));

        // stop when we hit the next contig's reads
        if(ctg.id != exl.match_chromosome()) break;

        const export_read_parse erp(exr);

        alignment al;
        al.pos=(exl.match_position()-1);
        al.is_fwd_strand=erp.is_fwd_strand;

        {   // read apath will be derived from the contig apath, but
            // first assert that there are no indels in the read's
            // alignment to the contig:
            const char* const md(exl.match_descriptor());
            export_md_to_apath(md,al.is_fwd_strand,al.path);
            assert((al.path.size()==1) && (al.path[0].type == ALIGNPATH::MATCH));
        }

        const bool is_ref_align(map_grouper_contig_read_to_genome(ctg,al));

        // this read could not be reference aligned (ie. it is entirely aligned within an insert)
        //
        if(! is_ref_align) continue;

        // TODO:rescue these reads as soft-clips if possible:
        //
        if(al.pos < 0) {
            std::string key;
            get_read_key_from_export_line(exl,key);
            log_os << "WARNING: skipping contig alignment for read: " << key << " due to negative alignment position: " << al.pos << "\n";
            continue;
        }

        // double-check:
        if(ALIGNPATH::is_apath_invalid(al.path,erp.read_size)) {
            log_os << "ERROR: Invalid contig read reference alignment: " << al;
            log_os << "\tread_size:" << erp.read_size << "\n";
            exr.report_state(log_os);
            exit(EXIT_FAILURE);
        }

        static const MAPLEVEL::index_t maplev(MAPLEVEL::UNKNOWN_MAPPED); // note this value is currently ignored for contig reads
        static const READ_ALIGN::index_t rat(READ_ALIGN::CONTIG);
        try {
            add_export_read_to_buffer(sppr,exl,erp.fwd_strand_read,erp.read_size,al,rat,
                                      ctg.chrom.c_str(),maplev,tmp_key_br,
                                      sample_no,contig_no,&contig_indels);
        } catch (...) {
            log_os << "\nException caught in add_export_read_to_buffer() while processing contig read alignment record:\n";
            exr.report_state(log_os);
            throw;
        }
    } while (exr.next());
}



void
process_contig(const starling_options& client_opt,
               const reference_contig_segment& ref,
               const grouper_contig& ctg,
               starling_pos_processor_base& sppr,
               const unsigned sample_no,
               const char* sample_label) {

    contig_no++;

#if 1//def DEBUG_STARLING
    log_os << "contig_no: " << contig_no << "\n";
    log_os << "contig: " << ctg;
#endif

    if(! is_valid_seq(ctg.seq.c_str())) {
        log_os << "ERROR: invalid sequence in contig: " << ctg << "\n";
        exit(EXIT_FAILURE);
    }

    static const INDEL_ALIGN_TYPE::index_t iat(INDEL_ALIGN_TYPE::CONTIG);
    const string_bam_seq bseq(ctg.seq);
    try {
        static const std::pair<bool,bool> edge_pin(std::make_pair(false,false));
        add_alignment_indels_to_sppr(client_opt.max_indel_size,ref,ctg,bseq,sppr,iat,contig_no,sample_no, edge_pin);
    } catch (...) {
        log_os << "\nException caught in add_alignment_indels_to_sppr() while processing contig";
        if(NULL != sample_label) log_os << " from sample '" << sample_label << "'";
        log_os << ":\n" << ctg;
        throw;
    }
}
