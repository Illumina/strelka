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

#include "somatic_result_set.hh"
#include "strelka_shared.hh"
#include "strelka_digt_states.hh"

#include "starling_common/starling_indel_call_pprob_digt.hh"

#include "boost/utility.hpp"


/// Object to manage the somatic indel probability model
///
/// This is formatted as an object only to help with pre-computation of the
/// prior.
///
struct somatic_indel_caller_grid : private boost::noncopyable
{
    explicit somatic_indel_caller_grid(const strelka_options& opt);

    void
    get_somatic_indel(
        const strelka_options& opt,
        const strelka_deriv_options& dopt,
        const starling_sample_options& normal_opt,
        const starling_sample_options& tumor_opt,
        const IndelKey& indelKey,
        const IndelData& indelData,
        const unsigned normalSampleIndex,
        const unsigned tumorSampleIndex,
        const bool is_use_alt_indel,
        somatic_indel_call& sindel) const;

private:
    blt_float_t _ln_som_match;
    blt_float_t _ln_som_mismatch;

    /// Germline genotype prior used for somatic indel calling
    blt_float_t _germlineGenotypeLogPrior[SOMATIC_DIGT::SIZE];
};
