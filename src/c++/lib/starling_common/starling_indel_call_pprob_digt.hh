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


/// get expected allele ratio for long indels
void
get_het_observed_allele_ratio(
    const unsigned read_length,
    const unsigned min_overlap,
    const IndelKey& indelKey,
    const double het_allele_ratio,
    double& log_ref_prob,
    double& log_indel_prob);


/// precalculates prior distributions based on theta value:
///
struct indel_digt_caller : private boost::noncopyable
{
    explicit
    indel_digt_caller(const double theta);

#if 0
    /// \brief call an indel @ pos by calculating the posterior probability
    /// of all possible genotypes for a diploid individual.
    ///
    void
    starling_indel_call_pprob_digt(
        const starling_base_options& client_opt,
        const starling_base_deriv_options& client_dopt,
        const starling_sample_options& sample_opt,
        const IndelKey& indelKey,
        const IndelSampleData& indelSampleData,
        const bool is_use_alt_indel,
        starling_diploid_indel& dindel) const;
#endif

    const double*
    lnprior_genomic(const bool is_haploid = false) const
    {
        return get_prior(is_haploid).genome;
    }

    // this prior isn't current used for single-sample indel calling
    // itself, but is available for indel_digt_caller clients:
    //
    const double*
    lnprior_polymorphic(const bool is_haploid = false) const
    {
        return get_prior(is_haploid).poly;
    }

    static
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

    static
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

    struct prior_group
    {
        double genome[STAR_DIINDEL::SIZE];
        double poly[STAR_DIINDEL::SIZE];
    };

private:

    const prior_group&
    get_prior(const bool is_haploid) const
    {
        return (is_haploid ? _lnprior_haploid : _lnprior);
    }

    prior_group _lnprior;
    prior_group _lnprior_haploid;
};
