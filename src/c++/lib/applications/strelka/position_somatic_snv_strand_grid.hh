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
///
/// variation on the original strawman snv caller -- implements a
/// compile-time specified grid in allele frequency space and requires
/// similar frequency as definition of non-somatic.
///
/// \author Chris Saunders
///

#pragma once

#include "somatic_result_set.hh"
#include "strelka_shared.hh"
#include "strelka_digt_states.hh"

#include "blt_common/position_snp_call_pprob_digt.hh"


/// Object used to pre-compute somatic snv priors
struct somatic_snv_caller_strand_grid
{
    explicit somatic_snv_caller_strand_grid(
        const strelka_options& opt);

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
    blt_float_t _contam_tolerance;
    blt_float_t _ln_csse_rate;
    blt_float_t _ln_sse_rate;

    blt_float_t _ln_som_match;
    blt_float_t _ln_som_mismatch;

    /// Germline genotype prior used for somatic SNV calling
    blt_float_t _germlineGenotypeLogPrior[SOMATIC_DIGT::SIZE];
};
