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
/// \author Chris Saunders, Sangtae Kim
///

#pragma once

#include "somatic_result_set.hh"
#include "blt_util/blt_types.hh"


/// \brief Calculates diploid germline variant prior probabilities used in the somatic variant model
///
/// Prior distribution is returned in log-space
void
calculateGermlineGenotypeLogPrior(
    const double theta,
    blt_float_t* germlineGenotypeLogPrior);


/// \brief Calculates somatic variant quality scores
///
/// \param[in] logSharedErrorRate Log of the rate at which the same (systematic) error is observed in the tumor and normal samples
/// \param[in] logSharedErrorRateComplement Log of 1-SharedErrorRate
void
calculate_result_set_grid(
    const blt_float_t ssnv_freq_ratio,
    const blt_float_t logSharedErrorRate,
    const blt_float_t logSharedErrorRateComplement,
    const blt_float_t* normal_lhood,
    const blt_float_t* tumor_lhood,
    const blt_float_t* germlineGenotypeLogPrior,
    const blt_float_t lnmatch,
    const blt_float_t lnmismatch,
    result_set& rs);
