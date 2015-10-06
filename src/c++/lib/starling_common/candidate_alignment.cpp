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

#include "candidate_alignment.hh"

#include "blt_util/align_path_util.hh"
#include "starling_common/starling_base_shared.hh"

#include <cassert>

#include <iostream>



std::ostream&
operator<<(std::ostream& os,
           const candidate_alignment& cal)
{
    os << "CANDIDATE_ALIGNMENT: " << cal.al;
    os << "CANDIDATE_ALIGNMENT_LEAD: " << cal.leading_indel_key;
    os << "CANDIDATE_ALIGNMENT_TRAIL: " << cal.trailing_indel_key;

    return os;
}



// get the keys of the indels present in the candidate alignment
//
void
get_alignment_indels(
    const candidate_alignment& cal,
    const unsigned max_indel_size,
    indel_set_t& indels)
{
    using namespace ALIGNPATH;

    indels.clear();

    const path_t& path(cal.al.path);
    unsigned path_index(0);
    unsigned read_offset(0);
    pos_t ref_head_pos(cal.al.pos);

    const std::pair<unsigned,unsigned> ends(get_match_edge_segments(cal.al.path));
    const unsigned aps(path.size());
    while (path_index<aps)
    {

        const bool is_edge_segment((path_index<ends.first) || (path_index>ends.second));
        const bool is_swap_start(is_segment_swap_start(path,path_index));

        const path_segment& ps(path[path_index]);

        assert(! (ps.type == SKIP));

        unsigned n_seg(1); // number of path_segments consumed
        if (is_swap_start)
        {
            const swap_info sinfo(path,path_index);
            n_seg=sinfo.n_seg;
        }

        if       (is_edge_segment)
        {
            // ignore all edge segments except INSERT/DELETE (this includes SWAP):
            if ((DELETE == ps.type) ||
                (INSERT == ps.type))
            {
                if (path_index<ends.first)
                {
                    assert(cal.leading_indel_key.type != INDEL::NONE);
                    indels.insert(cal.leading_indel_key);
                }
                else
                {
                    assert(cal.trailing_indel_key.type != INDEL::NONE);
                    indels.insert(cal.trailing_indel_key);
                }
            }

        }
        else if (is_swap_start)
        {
            const swap_info sinfo(path,path_index);

            const unsigned swap_size(std::max(sinfo.insert_length,sinfo.delete_length));
            if (swap_size <= max_indel_size)
            {
                indels.insert(indel_key(ref_head_pos,INDEL::SWAP,sinfo.insert_length,sinfo.delete_length));
            }
            else
            {
                indels.insert(indel_key(ref_head_pos,INDEL::BP_LEFT,swap_size));
                const pos_t right_pos(ref_head_pos+sinfo.delete_length);
                indels.insert(indel_key(right_pos,INDEL::BP_RIGHT,swap_size));
            }

        }
        else if (is_segment_type_indel(path[path_index].type))
        {
            if (ps.length <= max_indel_size)
            {
                const INDEL::index_t id( (ps.type==INSERT) ? INDEL::INSERT : INDEL::DELETE );
                indels.insert(indel_key(ref_head_pos,id,ps.length));
            }
            else
            {
                indels.insert(indel_key(ref_head_pos,INDEL::BP_LEFT,ps.length));
                const pos_t right_pos(ref_head_pos+((ps.type==INSERT) ? 0 : ps.length));
                indels.insert(indel_key(right_pos,INDEL::BP_RIGHT,ps.length));
            }

        }

        for (unsigned i(0); i<n_seg; ++i)
        {
            increment_path(path,path_index,read_offset,ref_head_pos);
        }
    }
}
