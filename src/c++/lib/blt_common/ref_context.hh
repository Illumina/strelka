// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
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

#include "blt_util/blt_types.hh"
#include "blt_util/reference_contig_segment.hh"


/// get size of hpol from the current position, given that this
/// is the left-most position in the hpol, otherwise, return 1
unsigned
get_left_shifted_hpol_size(
    const pos_t pos,
    const reference_contig_segment& ref);


/// Get the length of the longest homopolymer containing the current
/// position if this position can be treated as any base.
///
unsigned
get_snp_hpol_size(const pos_t pos,
                  const reference_contig_segment& ref);

/// find the largest homopolymer extending from pos where one
/// occurrence of an alternate base is allowed
///
unsigned
get_interrupted_hpol_size(
    const pos_t pos,
    const reference_contig_segment& ref);
