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

#include "alignment.hh"
#include "candidate_alignment.hh"


typedef std::vector<const candidate_alignment*> cal_pool_t;


/// given a set of high scoring alignments, compare them to identify ambiguous
/// segments on the alignment ends. Take the alignment for best_cal_id and
/// soft-clip any such ambiguous regions
void
get_clipped_alignment_from_cal_pool(
    const cal_pool_t& max_cal_pool,
    const unsigned best_cal_id,
    alignment& al);
