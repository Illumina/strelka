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
    const bool is_het_bias,
    const blt_float_t het_bias,
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

/// added by Sangtae to support 3 state model
void
get_diploid_gt_lhood_cached_simple(
    const snp_pos_info& pi,
    const unsigned ref_gt,
    blt_float_t* const lhood);

/// get lhood for nonstandard diploid het allele ratios
///
/// lhood ordering follows (undocumented) strelka conventions
///
/// \param hetresolution how many intermediates between 0-0.5 should we sample per half-axis?
void
get_diploid_het_grid_lhood_cached(
    const snp_pos_info& pi,
    const unsigned ref_gt,
    const unsigned hetResolution,
    blt_float_t* const lhood);
