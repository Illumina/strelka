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

private:
    std::vector<blt_float_t> _ln_freq_given_somatic;

    blt_float_t _ln_som_match;
    blt_float_t _ln_som_mismatch;
    const blt_float_t* _bare_lnprior;

    const std::vector<blt_float_t>&
    ln_freq_given_somatic() const
    {
        return _ln_freq_given_somatic;
    }

};
