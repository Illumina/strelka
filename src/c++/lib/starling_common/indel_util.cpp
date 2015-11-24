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

/// \file
///
/// \author Chris Saunders
///

#include "starling_common/indel_util.hh"



bool
is_indel_conflict(const indel_key& ik1,
                  const indel_key& ik2)
{

    // add one to the end_pos of all indels to prevent immediately
    // adjacent indels in the final alignments:
    pos_range pr1(ik1.open_pos_range());
    pr1.end_pos++;
    pos_range pr2(ik2.open_pos_range());
    pr2.end_pos++;

    return pr1.is_range_intersect(pr2);
}



bool
is_range_intersect_indel_breakpoints(const known_pos_range read_pr,
                                     const indel_key& ik)
{

    if (read_pr.is_range_intersect(pos_range(ik.pos,ik.pos))) return true;
    const pos_t rpos(ik.right_pos());
    if (ik.pos==rpos) return false;
    return (read_pr.is_range_intersect(pos_range(rpos,rpos)));
}



bool
is_range_adjacent_indel_breakpoints(const known_pos_range read_pr,
                                    const indel_key& ik)
{

    if (read_pr.is_range_intersect(pos_range(ik.pos-1,ik.pos+1))) return true;
    const pos_t rpos(ik.right_pos());
    if (ik.pos==rpos) return false;
    return (read_pr.is_range_intersect(pos_range(rpos-1,rpos+1)));
}
