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
/// \author Sangtae Kim
///

#include <blt_util/align_path_util.hh>
#include "ComplexAlleleUtil.hh"
#include "starling_read_util.hh"
#include "starling_pos_processor_base.hh"

unsigned
addComplexAlleleToSppr(
        const unsigned max_indel_size,
        const reference_contig_segment &ref,
        const alignment &al,
        const bam_seq_base &read_seq,
        starling_pos_processor_base& sppr,
        const INDEL_ALIGN_TYPE::index_t iat,
        const align_id_t id,
        const unsigned /*sample_no*/,
        const std::pair<bool, bool> &edge_pin,
        const bool is_mapq_zero)
{
    using namespace ALIGNPATH;

    const unsigned seq_len(read_seq.size());
//    log_os << rs << "\n";
    if (is_apath_invalid(al.path,seq_len))
    {
        std::ostringstream oss;
        oss << "ERROR: Can't handle alignment path '" << apath_to_cigar(al.path) << "' -- " << get_apath_invalid_reason(al.path,seq_len) << "\n";
        throw blt_exception(oss.str().c_str());
    }

    if (is_apath_starling_invalid(al.path))
    {
        std::ostringstream oss;
        oss << "ERROR: can't handle alignment path '" << apath_to_cigar(al.path) << "'\n";
        throw blt_exception(oss.str().c_str());
    }

    const rc_segment_bam_seq ref_bseq(ref);

    const std::pair<unsigned,unsigned> ends(get_match_edge_segments(al.path));

    pos_range valid_pr;
    get_valid_alignment_range(al,ref_bseq,read_seq,valid_pr);

    unsigned path_index(0);
    unsigned read_offset(0);
    pos_t ref_head_pos(al.pos);

    unsigned total_indel_ref_span_per_read(0);
    const unsigned aps(al.path.size());

    auto& active_region_detector(sppr.get_active_region_detector());

    while (path_index<aps)
    {
        const path_segment& ps(al.path[path_index]);
        const bool is_begin_edge(path_index<ends.first);
        const bool is_end_edge(path_index>ends.second);
        const bool is_edge_segment(is_begin_edge || is_end_edge);

        const bool is_swap_start(is_segment_swap_start(al.path,path_index));

        assert(ps.type != SKIP);

        // fix STARKA-43:
        // assert(! (is_edge_segment && is_swap_start));

        IndelObservation obs;
        obs.data.iat = iat;
        obs.data.id = id;

        if (! is_segment_align_match(ps.type))
        {
//            log_os << al.path <<" \n";
            pos_range indel_read_pr;
            indel_read_pr.set_begin_pos((read_offset==0) ? 0 : (read_offset-1));

            unsigned rlen(0);
            if       (is_swap_start)
            {
                const swap_info sinfo(al.path,path_index);
                rlen=sinfo.insert_length;

                if (sinfo.delete_length<=max_indel_size)
                {
                    total_indel_ref_span_per_read += sinfo.delete_length;
                }
            }
            else if (is_segment_type_read_length(ps.type))
            {
                rlen=ps.length;
            }
            else
            {
                if (ps.type == DELETE)
                {
                    if (ps.length <= max_indel_size)
                    {
                        total_indel_ref_span_per_read += ps.length;
                    }
                }
            }
            indel_read_pr.set_end_pos(std::min(seq_len,read_offset+1+rlen));
            if (! valid_pr.is_superset_of(indel_read_pr)) obs.data.is_noise=true;
        }

        unsigned n_seg(1); // number of path segments consumed
        if (is_edge_segment)
        {
            // is this indel occurring on a pinned edge (ie against an exon?)
            const bool is_pinned_indel((is_begin_edge && edge_pin.first) ||
                                       (is_end_edge && edge_pin.second));

            // edge inserts are allowed for intron adjacent and grouper reads, edge deletions for intron adjacent only:
            if (ps.type == INSERT)
            {
//                process_edge_insert(max_indel_size,al.path,read_seq,
//                                    sppr,obs,sample_no,
//                                    path_index,read_offset,ref_head_pos,
//                                    is_pinned_indel);
            }
            else if (ps.type == DELETE)
            {
                if (is_pinned_indel)
                {
//                    process_edge_delete(max_indel_size,al.path,read_seq,
//                                        sppr,obs,sample_no,
//                                        path_index,read_offset,ref_head_pos,
//                                        is_pinned_indel);
                }
            }
        }
        else if (is_swap_start)
        {
//            n_seg = process_swap(max_indel_size,al.path,read_seq,
//                                 sppr,obs,sample_no,
//                                 path_index,read_offset,ref_head_pos);

        }
        else if (is_segment_type_indel(al.path[path_index].type))
        {
//            process_simple_indel(max_indel_size,al.path,read_seq,
//                                 sppr,obs,sample_no,
//                                 path_index,read_offset,ref_head_pos);
        }
        else if (!is_mapq_zero && sppr.is_active_region_detector_enabled() && is_segment_align_match(ps.type))
        {
            // detect active regions (match/mismatch)
            for (unsigned j(0); j < ps.length; ++j)
            {
                const pos_t ref_pos(ref_head_pos + static_cast<pos_t>(j));
                pos_t read_pos = read_offset + j;

                char base_char = read_seq.get_char(read_pos);

                if (ref.get_base(ref_pos) != base_char)
                {
                    active_region_detector.insertMismatch(id, ref_pos, base_char);
                }
                else
                {
                    active_region_detector.insertMatch(id, ref_pos, base_char);
                }
            }
        }

        if (sppr.is_active_region_detector_enabled() && !is_mapq_zero && obs.key.type != INDEL::NONE)
        {
            // detect active region (indel)
            active_region_detector.insertIndel(obs);
        }

        for (unsigned i(0); i<n_seg; ++i)
        {
            increment_path(al.path,path_index,read_offset,ref_head_pos);
        }
    }

    return total_indel_ref_span_per_read;
}

