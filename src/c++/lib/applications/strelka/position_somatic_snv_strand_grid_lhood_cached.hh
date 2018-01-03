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

/// \author Chris Saunders
/// \author Sangtae Kim
///

#pragma once

#include "blt_common/blt_shared.hh"
#include "blt_common/snp_pos_info.hh"
#include "strelka_common/het_ratio_cache.hh"

void
get_diploid_gt_lhood_cached_simple(
    const snp_pos_info& pi,
    const unsigned ref_gt,
    blt_float_t* const lhood);

void
get_diploid_het_grid_lhood_cached(
    const snp_pos_info& pi,
    const unsigned ref_gt,
    const unsigned hetResolution,
    blt_float_t* const lhood);

void
get_strand_ratio_lhood_spi(
    const snp_pos_info& pi,
    const unsigned ref_gt,
    const blt_float_t het_ratio,
    const unsigned het_ratio_index,
    het_ratio_cache<2>& hrcache,
    blt_float_t* lhood);
