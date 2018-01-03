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

#pragma once

#include "blt_util/pos_range.hh"
#include "starling_common/indel.hh"


/// Return true if the two indels are overlapping or adjacent
bool
is_indel_conflict(
    const IndelKey& indelKey1,
    const IndelKey& indelKey2);


/// check if a range (representing an alignment, etc) intersects the breakends of a
/// candidate indel.
///
/// note: (1) intersection does not include adjacency
///       (2) this is breakpoint intersection, so for example a spanning deletion will not be included
///
bool
is_range_intersect_indel_breakpoints(
    const known_pos_range read_pr,
    const IndelKey& indelKey);


/// check if a range (representing an alignment, etc) intersects or is
/// adjacent to the breakends of a candidate indel
///
bool
is_range_adjacent_indel_breakpoints(
    const known_pos_range read_pr,
    const IndelKey& indelKey);
