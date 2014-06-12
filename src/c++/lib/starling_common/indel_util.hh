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

/// \file
///
/// \author Chris Saunders
///

#pragma once

#include "blt_util/pos_range.hh"
#include "starling_common/indel.hh"


// do these two indels overlap?:
bool
is_indel_conflict(const indel_key& ik1,
                  const indel_key& ik2);


/// check if a range (representing an alignment, etc) intersects the breakends of a
/// candidate indel.
///
/// note: (1) intersection does not include adjacency
///       (2) this is breakpoint intersection, so for example a spanning deletion will not be included
///
bool
is_range_intersect_indel_breakpoints(const known_pos_range read_pr,
                                     const indel_key& ik);


/// check if a range (representing an alignment, etc) intersects or is
/// adjacent to the breakends of a candidate indel
///
bool
is_range_adjacent_indel_breakpoints(const known_pos_range read_pr,
                                    const indel_key& ik);
