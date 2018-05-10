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

#include "starling_common/starling_base_shared.hh"
#include "starling_common/starling_diploid_indel.hh"
#include "starling_common/AlleleReportInfo.hh"
#include "starling_common/indel.hh"

#include <iosfwd>


/// total the path likelihoods of ref,indel and alt_indel states
///
/// utility used for indel filtration...
///
void
get_sum_path_pprob(
    const starling_base_deriv_options& dopt,
    const IndelSampleData& indelSampleData,
    const bool is_include_tier2,
    const bool is_use_alt_indel,
    ReadPathScores& total_pprob,
    const bool is_init_total = true);


/// get expected allele ratio between two indel alleles (where one allele must be reference)
///
/// this  only serves an important role for indels which are long relative to the read size
/// For sites and single breakpoints this is expected to match the sample allele ratio.
///
/// Currently assuming second allele is the reference, generalization TBD
///
void
get_het_observed_allele_ratio(
    const unsigned read_length,
    const unsigned min_overlap,
    const IndelKey& indelKey,
    const double het_allele_ratio,
    double& log_ref_prob,
    double& log_indel_prob);

/// get heterozygous indel likelihoods for non-canonical allele ratios
void
get_high_low_het_ratio_lhood(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const starling_sample_options& sample_opt,
    const IndelKey& indelKey,
    const IndelSampleData& indelSampleData,
    const double het_ratio,
    const bool is_include_tier2,
    const bool is_use_alt_indel,
    double& het_lhood_high,
    double& het_lhood_low);

/// get diploid indel likelihoods
///
/// this is the old function for diploid likelihoods, which assumes we're evaluating one indel at a time without
/// considering overlap.
void
get_indel_digt_lhood(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const starling_sample_options& sample_opt,
    const IndelKey& indelKey,
    const IndelSampleData& indelSampleData,
    const bool is_include_tier2,
    const bool is_use_alt_indel,
    double* const lhood);
