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

#pragma once


//
// Breakpoints refer to both insertions and deletions which exceed
// varling's MAX_INDEL_SIZE. In this case we try to account for the
// breakpoint during snp-calling, but do not attempt to call the indel
// using the same methods used for small indels.
//
// Positions: Internally all positions are stored using zero-indexed
// position numbers.  For small indels and left breakpoints, we store
// the position of the first affected position. For right breakpoints,
// we store the first positions *after* the last affected
// position. Positions are stored in this manner so that the indels
// follow the starling range convention
//
// Large indel breakpoints are expected to be rare and to possibly
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
// The Swap type is initially used to represent combined
// insertion/deletion events, but will be used in the future for
// general alternate haplotypes.
//
// The "NONE" type is used for some indel lookup methods, because it sorts
// ahead of all other types at a given position:
//
namespace INDEL {
enum index_t {
    NONE,
    INSERT,
    DELETE,
    BP_LEFT,
    BP_RIGHT,
    SWAP
};

inline
const char*
get_index_label(index_t id) {
    switch(id) {
    case NONE:        return "NONE";
    case INSERT:      return "INSERT";
    case DELETE:      return "DELETE";
    case BP_LEFT:  return "BP_LEFT";
    case BP_RIGHT: return "BP_RIGHT";
    case SWAP:        return "SWAP";
    default:          return "UNKNOWN";
    }
}
}

