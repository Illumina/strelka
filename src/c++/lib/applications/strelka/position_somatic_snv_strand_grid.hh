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

/// variation on the original strawman snv caller -- implements a
/// compile-time specified grid in allele frequency space and requires
/// similar frequency as definition of non-somatic.
///

/// \author Chris Saunders
///

#pragma once

#include "extended_pos_data.hh"
#include "position_somatic_snv_grid_shared.hh"
#include "position_somatic_snv_strand_grid_states.hh"
#include "strelka_shared.hh"

#include "blt_common/position_snp_call_pprob_digt.hh"


// object used to pre-compute priors:
struct somatic_snv_caller_strand_grid
{
    somatic_snv_caller_strand_grid(
        const strelka_options& opt,
        const pprob_digt_caller& pd_caller);

    //
    void
    position_somatic_snv_call(
        const extended_pos_info& normal_epi,
        const extended_pos_info& tumor_epi,
        const extended_pos_info* normal_epi_t2_ptr,
        const extended_pos_info* tumor_epi_t2_ptr,
        const bool isComputeNonSomatic,
        somatic_snv_genotype_grid& sgt) const;

    // compute a lot of prior information for various alternate
    // versions of the method -- we don't actually need all of this for any one computation:
    //
    struct prior_set
    {
        prior_set()
            : normal(DIGT_SGRID::SIZE)
            , somatic_marginal(DIGT_SGRID::SIZE)
            , normal_poly(DIGT_SGRID::SIZE)
            , somatic_marginal_poly(DIGT_SGRID::SIZE)
            , normal_nostrand(DIGT_SGRID::SIZE)
            , normal_poly_nostrand(DIGT_SGRID::SIZE)
        {}

        typedef std::vector<blt_float_t> prior_t;

        prior_t normal;
        prior_t somatic_marginal;
        prior_t normal_poly;
        prior_t somatic_marginal_poly;

        // added to support somatic gVCF:
        prior_t normal_nostrand;
        prior_t normal_poly_nostrand;
    };

private:

    const prior_set&
    get_prior_set(const unsigned ref_id) const
    {
        return _lnprior[ref_id];
    }

    const std::vector<blt_float_t>&
    lnprior_genomic(const unsigned ref_id) const
    {
        return _lnprior[ref_id].normal;
    }

    const std::vector<blt_float_t>&
    lnprior_polymorphic(const unsigned ref_id) const
    {
        return _lnprior[ref_id].normal_poly;
    }

    const strelka_options& _opt;
    std::array<prior_set,N_BASE+1> _lnprior;
    blt_float_t _ln_som_match;
    blt_float_t _ln_som_mismatch;
};
