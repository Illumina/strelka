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

///
/// \author Chris Saunders
///

#pragma once


// Breakpoints refer to both insertions and deletions which exceed
// the method's maxIndelSize. In this case we try to account for the
// breakpoint during realignment (to improve snp-calling), but do not
// attempt to call the indel using the same methods used for small indels.
//
// Positions: Internally all positions are stored using zero-indexed
// position numbers.  For small indels and left breakpoints, we store
// the position of the first affected position. For right breakpoints,
// we store the first positions *after* the last affected
// position. Positions are stored in this manner so that the indels
// follow Strelka's range convention
//
// Large indels are expected to be rare and to possibly
// have an unknown size (especially for insertions). Thus the indel
// length is not used in this case. To accommodate multiple different
// lengths at a single location, an advisory length may be inserted
// but it will not be used as part of any alignment calculation. For a
// breakpoint corresponding to a large deletion, "seq" is expected to
// represent the sequence on the other side of the breakpoint,
// starling will not look it up. In this way the break-point for *any*
// type of large-scale event can be accommodated for snp-calling. The
// seq convention for a breakpoint is to set the seq size equal to
// MAX_INDEL_SIZE
//
// Skips (introns) are *not* added here. They are treated as invariant
// alignment elements in the original read.
//
// The "NONE" type is used for some indel lookup methods, because it sorts
// ahead of all other types at a given position:
//
namespace INDEL
{
enum index_t
{
    NONE,
    INDEL,
    MISMATCH,
    BP_LEFT,
    BP_RIGHT,
};

inline
const char*
get_index_label(index_t id)
{
    switch (id)
    {
    case NONE:
        return "NONE";
    case INDEL:
        return "INDEL";
    case MISMATCH:
        return "MISMATCH";
    case BP_LEFT:
        return "BP_LEFT";
    case BP_RIGHT:
        return "BP_RIGHT";
    default:
        return "UNKNOWN";
    }
}
}

