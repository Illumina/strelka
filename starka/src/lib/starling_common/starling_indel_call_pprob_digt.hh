// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
//

/// \file

/// \author Chris Saunders
///
#ifndef __STARLING_INDEL_CALL_PPROB_DIGT_HH
#define __STARLING_INDEL_CALL_PPROB_DIGT_HH

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


// debuging output -- produces labels
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
struct indel_digt_caller : private boost::noncopyable {
    
    indel_digt_caller(const double theta);

    // \brief call an indel @ pos by calculating the posterior probability
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
    lnprior_genomic() const { return _lnprior_genomic; }

    // this prior isn't current used for single-sample indel calling
    // itself, but is available for indel_digt_caller clients:
    //
    const double*
    lnprior_polymorphic() const { return _lnprior_polymorphic; }

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

private:
    double _lnprior_genomic[STAR_DIINDEL::SIZE];
    double _lnprior_polymorphic[STAR_DIINDEL::SIZE];
};


#endif
