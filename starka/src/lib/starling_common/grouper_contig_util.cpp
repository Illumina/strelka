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


#include "blt_util/blt_exception.hh"
#include "blt_util/export_stream_reader.hh"
#include "blt_util/io_util.hh"
#include "blt_util/log.hh"
#include "starling_common/align_path_util.hh"
#include "starling_common/grouper_contig_util.hh"
#include "starling_common/starling_shared.hh"

#include "boost/lexical_cast.hpp"

#include <cassert>
#include <cstdio> // for EOF...

#include <fstream>
#include <iostream>
#include <sstream>



// reads which extend off the end of the contig on either side are
// soft-clipped
//
bool
map_grouper_contig_read_to_genome(const grouper_contig& ctg,
                                  alignment& read_al) {

    using namespace ALIGNPATH;

    // Read alignment is expected to perfectly align to the contig
    // (but soft-clip is allowed):
    //
    const unsigned rp_size(read_al.path.size());
    assert(rp_size>=1);
    for (unsigned i(0); i<rp_size; ++i) {
        const align_t ali(read_al.path[i].type);
        assert((ali==MATCH) || (ali==SOFT_CLIP));
    }

    const unsigned al_ref_length(apath_ref_length(read_al.path));

    pos_t read_begin_pos(read_al.pos);
    const pos_t read_end_pos(read_begin_pos+al_ref_length);

    assert(read_end_pos>=0);  // else read doesn't align within the contig

    path_t new_path;

    // add leading soft-clip segment if it exists:
    if (read_al.path[0].type == SOFT_CLIP) {
        new_path.push_back(read_al.path[0].type);
    }

    // mark any sequence hanging off the end of the contig as soft-clipped:
    if (read_begin_pos<0) {
        path_segment seg;
        seg.type = SOFT_CLIP;
        seg.length = -read_begin_pos;
        new_path.push_back(seg);
        read_begin_pos=0;
    }

    bool is_read_path_begin(false);
    unsigned contig_segment_offset(0);
    pos_t ref_segment_end_pos(ctg.pos);
    const unsigned cpath_size(ctg.path.size());

    unsigned path_index(0);
    while (path_index<cpath_size) {
        path_segment seg(ctg.path[path_index]);

        const bool is_swap_start(is_segment_swap_start(ctg.path,path_index));

        const pos_t contig_last_segment_offset(contig_segment_offset);
        const pos_t ref_last_segment_end_pos(ref_segment_end_pos);

        unsigned n_seg(1);
        std::auto_ptr<swap_info> sinfo_ptr;
        if (is_swap_start) {
            sinfo_ptr.reset(new swap_info(ctg.path,path_index));
            n_seg=sinfo_ptr->n_seg;
        }

        for (unsigned i(0); i<n_seg; ++i) {
            increment_path(ctg.path,path_index,
                           contig_segment_offset,
                           ref_segment_end_pos);
        }

        bool is_edge_segment(false);
        unsigned read_segment_reduction(0);

        if (! is_read_path_begin) {
            if (read_begin_pos >= static_cast<pos_t>(contig_segment_offset)) continue;

            is_read_path_begin=true;
            is_edge_segment=true;

            // read_begin_pos is in contig coordinates, so "shift"
            // refers to the additional offset to start the read from
            // the begining of the last contig CIGAR segment
            //
            const pos_t shift(read_begin_pos-contig_last_segment_offset);

            if (is_swap_start) {
                read_al.pos = ref_segment_end_pos;
                assert(shift < static_cast<pos_t>(read_segment_reduction+sinfo_ptr->insert_length));
                read_segment_reduction += shift;
            } else {
                read_al.pos = std::min(ref_last_segment_end_pos+shift,ref_segment_end_pos);
                assert(shift < static_cast<pos_t>(seg.length));
                seg.length -= shift;
            }
        }

        const bool is_final_seg(read_end_pos <= static_cast<pos_t>(contig_segment_offset));
        if (is_final_seg) {
            is_edge_segment=true;
            const pos_t shift(contig_segment_offset-read_end_pos);
            if (is_swap_start) {
                assert(shift < static_cast<pos_t>(read_segment_reduction+sinfo_ptr->insert_length));
                read_segment_reduction += shift;
            } else {
                assert(shift < static_cast<pos_t>(seg.length));
                seg.length -= shift;
            }
        }

        if (is_swap_start) {
            assert(sinfo_ptr->insert_length > read_segment_reduction);
            new_path.push_back(path_segment(INSERT,(sinfo_ptr->insert_length-read_segment_reduction)));
            if (! is_edge_segment) {
                new_path.push_back(path_segment(DELETE,sinfo_ptr->delete_length));
            }
        } else {
            new_path.push_back(seg);
        }

        if (is_final_seg) break;
    }


    assert(is_read_path_begin); // else read doesn't align within the contig

    // add trailing soft-clip if it exists:
    if ((rp_size > 1) and (read_al.path[rp_size-1].type == SOFT_CLIP)) {
        new_path.push_back(read_al.path[rp_size-1]);
    }

    // mark any read segment hanging off the end of the contig as soft-clipped:
    const unsigned al_read_length(apath_read_length(read_al.path));
    const unsigned new_read_length(apath_read_length(new_path));
    if (new_read_length < al_read_length) {
        path_segment seg;
        seg.type = SOFT_CLIP;
        seg.length = (al_read_length-new_read_length);
        new_path.push_back(seg);
    }

    // copy and condense new path, check that a MATCH segment exists:
    bool is_genome_align(false);
    {
        path_t& rp(read_al.path);
        rp.clear();
        align_t last_type(NONE);
        const unsigned nps(new_path.size());
        for (unsigned i(0); i<nps; ++i) {
            const path_segment& ps(new_path[i]);
            if (ps.type == MATCH) is_genome_align=true;
            if (ps.type != last_type) {
                rp.push_back(ps);
                last_type=ps.type;
            } else {
                rp.back().length += ps.length;
            }
        }
    }
    return is_genome_align;
}



static
void
bad_header_error(const std::string& header) {
    std::ostringstream oss;
    oss << "ERROR: Unexpected format in GROUPER contig header: " << header << "\n";
    throw blt_exception(oss.str().c_str());
}



bool
get_next_contig(std::istream& is,
                grouper_contig& ctg) {

    ctg.clear();

    if (! is) return false;

    // get header:
    {
        const int c(is.peek());
        if (c==EOF) return false;
        if (c != '>') {
            throw blt_exception("ERROR: Unexpected format in GROUPER contig file\n");
        }

        std::string header;
        if (not std::getline(is,header)) {
            throw blt_exception("ERROR: Unexpected format in GROUPER contig file\n");
        }

        static const unsigned header_field_count(5);
        std::string::size_type begin_pos,pos(1);
        for (unsigned i(0); i<header_field_count; ++i) {
            static const char header_delim[] = "| \t\n\r";

            begin_pos=header.find_first_not_of(header_delim,pos);
            pos=header.find_first_of(header_delim,begin_pos);

            if (std::string::npos == begin_pos) bad_header_error(header);

            if       (i==0) {
                ctg.id=header.substr(begin_pos,pos-begin_pos);
            } else if (i==1) {
                ctg.chrom=header.substr(begin_pos,pos-begin_pos);
            } else if (i==2) {
                ctg.pos=boost::lexical_cast<pos_t>(header.substr(begin_pos,pos-begin_pos))-1;
            } else if (i==3) {
                const std::string cigar(header.substr(begin_pos,pos-begin_pos));
                cigar_to_apath(cigar.c_str(),ctg.path);
            } else if (i==4) {
                //intentional throwaway:
                //pos_t indel_begin_pos=boost::lexical_cast<pos_t>(header.substr(begin_pos,pos-begin_pos))-1;

            } else {
                bad_header_error(header);
            }
        }

        /// \todo TODO -- handle this case for circular chromosome support:
        assert(ctg.pos>=0);
    }

    // read seq:
    {
        if (is.peek() == '>') {
            throw blt_exception("ERROR: Unexpected format in GROUPER contig file.\n");
        }

        static const char seq_delim[] = " \t\n\r";
        std::string line_buff;

        std::string::size_type begin_pos,pos(1);
        std::string& seq(ctg.seq);
        while (std::getline(is,line_buff)) {
            begin_pos=line_buff.find_first_not_of(seq_delim);
            pos=line_buff.find_last_not_of(seq_delim)+1;

            if (ctg.is_usable) {
                seq += line_buff.substr(begin_pos,pos-begin_pos);
                if (seq.size() > MAX_CONTIG_SIZE) {
                    log_os << "WARNING: GROUPER contig: '" << ctg.id << "' exceeds maximum contig size. Skipping...\n";
                    ctg.is_usable=false;
                }
            }
            if (is.peek() == '>') break;
        }
        return true;
    }
}



contig_data_manager::
contig_data_manager(const std::string& contig_filename,
                    const std::string& contig_read_filename)
    : _contig_isp(new std::ifstream)
    , _contig_read_isp(new std::ifstream) {

    if (contig_filename.empty() != contig_read_filename.empty()) {
        log_os << "ERROR: contig and contig reads file must specified together or not at all\n";
        exit(EXIT_FAILURE);
    }

    std::ifstream& contig_is(*_contig_isp.get());
    std::ifstream& contig_read_is(*_contig_read_isp.get());

    if (not contig_filename.empty()) {
        open_ifstream(contig_is,contig_filename.c_str());
        open_ifstream(contig_read_is,contig_read_filename.c_str());
    }

    _contig_read_exrp.reset(new export_stream_reader(contig_read_is,contig_filename.c_str()));
    _creaderp.reset(new contig_reader(contig_is));
}



// leave this in cpp file so that auto_ptr links correctly:
contig_data_manager::
~contig_data_manager() {}
