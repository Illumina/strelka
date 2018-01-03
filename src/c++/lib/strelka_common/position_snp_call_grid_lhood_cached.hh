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
///

#pragma once

#include "het_ratio_cache.hh"

#include "blt_common/blt_shared.hh"
#include "blt_common/snp_pos_info.hh"


/// get standard diploid snp lhood's
///
/// note this is similar to starling's version, except by
/// excluding dependent_eprob values we can cache intermediate
/// values for each qscore and run this more quickly:
///
/// note lhood is expected to follow the standard genotype order defined
/// in blt_util/digt.hh
///
void
get_diploid_gt_lhood_cached(
    const blt_options& opt,
    const snp_pos_info& pi,
    const bool useHetVariantFrequencyExtension,
    const blt_float_t hetVariantFrequencyExtension,
    blt_float_t* const lhood);


/// same as above but with het_bias conveniently taken out of the interface
inline
void
get_diploid_gt_lhood_cached(
    const blt_options& opt,
    const snp_pos_info& pi,
    blt_float_t* const lhood)
{
    get_diploid_gt_lhood_cached(opt, pi, false, 0, lhood);
}


/// get lhood for nonstandard diploid het allele ratios
///
/// lhood ordering follows (undocumented) strelka conventions
///
/// \param hetresolution how many intermediates between 0-0.5 should we sample per half-axis?
void
get_diploid_het_grid_lhood_cached(
    const snp_pos_info& pi,
    const unsigned hetResolution,
    blt_float_t* const lhood);
