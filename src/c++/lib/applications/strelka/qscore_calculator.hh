// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
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
/// \author Chris Saunders, Sangtae Kim
///

#pragma once

#include "somatic_result_set.hh"
#include "strelka_shared.hh"
#include "strelka_digt_states.hh"

#include "blt_util/math_util.hh"
#include "blt_common/position_snp_call_pprob_digt.hh"

#include <cmath>


// calculates diploid prior probabilities
void
calculate_bare_lnprior(const double theta,
                       blt_float_t* bare_lnprior);

// calculates qscores
void
calculate_result_set_grid(
    const blt_float_t ssnv_freq_ratio,
    const blt_float_t ln_se_rate,   // ln (shared_error_rate)
    const blt_float_t ln_cse_rate,  // ln (1 - shared_error_rate)
    const blt_float_t* normal_lhood,
    const blt_float_t* tumor_lhood,
    const blt_float_t* bare_lnprior,
    const blt_float_t lnmatch,
    const blt_float_t lnmismatch,
    result_set& rs
);
