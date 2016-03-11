#pragma once

#include "position_somatic_snv_grid_shared.hh"
#include "position_somatic_snv_strand_grid_states.hh"
#include "strelka_shared.hh"

#include "blt_common/position_snp_call_pprob_digt.hh"

void
calculate_result_set_grid(
        const bool isComputeNonSomatic,
        const blt_float_t* normal_lhood,
        const blt_float_t* tumor_lhood,
        const std::vector<blt_float_t>& ln_freq_given_somatic,
        const blt_float_t* bare_lnprior_normal,
        const blt_float_t lnmatch,
        const blt_float_t lnmismatch,
        const bool is_forced_output,
        result_set& rs
        );
