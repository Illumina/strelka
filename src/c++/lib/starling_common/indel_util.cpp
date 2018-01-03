//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
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
/// \author Chris Saunders
///

#include "starling_common/indel_util.hh"



bool
is_indel_conflict(
    const IndelKey& indelKey1,
    const IndelKey& indelKey2)
{
    // if both indels are not mismatches
    // add one to the end_pos of all indels to prevent immediately
    // adjacent indels in the final alignments:

    bool isEitherMismatch(indelKey1.isMismatch() || indelKey2.isMismatch());
    pos_range pr1(indelKey1.open_pos_range());
    if (!isEitherMismatch) pr1.end_pos++;

    pos_range pr2(indelKey2.open_pos_range());
    if (!isEitherMismatch) pr2.end_pos++;

    return pr1.is_range_intersect(pr2);
}

/// return true if the range intersects the left or the right indel position
bool
is_range_intersect_indel_breakpoints(
    const known_pos_range read_pr,
    const IndelKey& indelKey)
{
    if (indelKey.isMismatch())
    {
        assert(indelKey.delete_length() == 1);
        return read_pr.is_pos_intersect(indelKey.pos);
    }

    if (read_pr.is_range_intersect(pos_range(indelKey.pos,indelKey.pos))) return true;
    const pos_t rpos(indelKey.right_pos());
    if (indelKey.pos==rpos) return false;
    return (read_pr.is_range_intersect(pos_range(rpos,rpos)));
}



bool
is_range_adjacent_indel_breakpoints(
    const known_pos_range read_pr,
    const IndelKey& indelKey)
{
    if (read_pr.is_range_intersect(pos_range(indelKey.pos-1,indelKey.pos+1))) return true;
    const pos_t rpos(indelKey.right_pos());
    if (indelKey.pos==rpos) return false;
    return (read_pr.is_range_intersect(pos_range(rpos-1,rpos+1)));
}
