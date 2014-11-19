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

///
/// \author Chris Saunders
///

#pragma once

#include "starling_common/starling_diploid_indel.hh"
#include "starling_common/starling_indel_report_info.hh"

#include "starling_common/indel.hh"
#include "starling_common/starling_shared.hh"

#include <iosfwd>


// total the path likelihoods of ref,indel and alt_indel states
//
// utility used for indel filtration...
//
void
get_sum_path_pprob(const starling_deriv_options& dopt,
                   const indel_data& id,
                   const bool is_include_tier2,
                   const bool is_use_alt_indel,
                   read_path_scores& total_pprob,
                   const bool is_init_total = true);


// Debugging output -- produces labels
//
void
write_starling_diploid_indel(const starling_diploid_indel& dgt,
                             const starling_indel_report_info& iri,
                             const starling_indel_sample_report_info& isri,
                             std::ostream& os);

// file output -- data only
//
void
write_starling_diploid_indel_file(const starling_diploid_indel& dgt,
                                  const starling_indel_report_info& iri,
                                  const starling_indel_sample_report_info& isri,
                                  std::ostream& os);


// Use caller object to precalculate prior distributions based on
// theta value:
//
struct indel_digt_caller : private boost::noncopyable
{
    explicit
    indel_digt_caller(const double theta);

    /// \brief call an indel @ pos by calculating the posterior probability
    /// of all possible genotypes for a diploid individual.
    ///
    void
    starling_indel_call_pprob_digt(const starling_options& client_opt,
                                   const starling_deriv_options& client_dopt,
                                   const starling_sample_options& sample_opt,
                                   const double indel_error_prob,
                                   const double ref_error_prob,
                                   const indel_key& ik,
                                   const indel_data& id,
                                   const bool is_use_alt_indel,
                                   starling_diploid_indel& dindel) const;

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
    get_high_low_het_ratio_lhood(const starling_options& opt,
                                 const starling_deriv_options& dopt,
                                 const starling_sample_options& sample_opt,
                                 const double indel_error_lnp,
                                 const double indel_real_lnp,
                                 const double ref_error_lnp,
                                 const double ref_real_lnp,
                                 const indel_key& ik,
                                 const indel_data& id,
                                 const double het_ratio,
                                 const bool is_include_tier2,
                                 const bool is_use_alt_indel,
                                 double& het_lhood_high,
                                 double& het_lhood_low);

    static
    void
    get_indel_digt_lhood(const starling_options& opt,
                         const starling_deriv_options& dopt,
                         const starling_sample_options& sample_opt,
                         const double indel_error_prob,
                         const double ref_error_prob,
                         const indel_key& ik,
                         const indel_data& id,
                         const bool is_het_bias,
                         const double het_bias,
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
