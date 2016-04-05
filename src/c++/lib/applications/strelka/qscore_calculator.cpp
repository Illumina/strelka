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
/// \author Sangtae Kim
///


#include "qscore_calculator.hh"
#include "somatic_call_shared.hh"
#include "blt_common/snp_util.hh"
#include "blt_util/log.hh"
#include "blt_util/math_util.hh"
#include "blt_util/prob_util.hh"
#include "blt_util/seq_util.hh"
#include "strelka_common/position_snp_call_grid_lhood_cached.hh"

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <map>

void
calculate_bare_lnprior(const double theta,
        blt_float_t *bare_lnprior)
{
    bare_lnprior[SOMATIC_DIGT::REF] = (blt_float_t) log1p_switch(-(3.*theta)/2.);
    bare_lnprior[SOMATIC_DIGT::HOM] = (blt_float_t) std::log(theta/2.);
    bare_lnprior[SOMATIC_DIGT::HET] = (blt_float_t) std::log(theta);
}

static const blt_float_t neg_inf = -std::numeric_limits<float>::infinity();
static const blt_float_t ln_one_half(std::log(1./2.));
static const blt_float_t log_error_mod = -std::log(static_cast<double>(DIGT_GRID::PRESTRAND_SIZE-1));

void
calculate_result_set_grid(
        const blt_float_t contam_tolerance,
        const blt_float_t ln_se_rate,   // ln (shared_error_rate)
        const blt_float_t ln_cse_rate,  // ln (1 - shared_error_rate)
        const blt_float_t* normal_lhood,
        const blt_float_t* tumor_lhood,
        const blt_float_t* bare_lnprior,
        const blt_float_t lnmatch,
        const blt_float_t lnmismatch,
        result_set& rs
        )
{
    double log_post_prob[SOMATIC_DIGT::SIZE][SOMATIC_STATE::SIZE];
    double max_log_prob = neg_inf;

    rs.max_gt = 0;

    for (unsigned ngt(0); ngt<SOMATIC_DIGT::SIZE; ++ngt)
    {
        for (unsigned tgt(0); tgt<SOMATIC_STATE::SIZE; ++tgt) // 0: non-somatic, 1: somatic
        {
            double max_log_sum = neg_inf;
            double log_sum[DDIGT_GRID::PRESTRAND_SIZE];
            int index = 0;

            for (unsigned tumor_freq_index(0); tumor_freq_index<DIGT_GRID::PRESTRAND_SIZE; ++tumor_freq_index)
            {
                blt_float_t tumor_freq = DIGT_GRID::get_fraction_from_index(tumor_freq_index);

                // this flag is set if tumor_freq is large.
                // if it is set, normal_freq==DIGT_GRID::RATIO_INCREMENT is considered as "canonical" frequency rather than noise
                bool consider_norm_contam = contam_tolerance*tumor_freq >= DIGT_GRID::RATIO_INCREMENT;

                for (unsigned normal_freq_index(0); normal_freq_index<DIGT_GRID::PRESTRAND_SIZE; ++normal_freq_index)
                {
                    double lprior_freq;

                    if (tgt == 0)    // non-somatic
                    {
                       if (normal_freq_index != tumor_freq_index) continue; // P(fn != ft | Gn = Gt) = 0

                       lprior_freq = (normal_freq_index == ngt) ? ln_cse_rate : ln_se_rate+log_error_mod;
                    }
                    else    // somatic
                    {
                       if (normal_freq_index == tumor_freq_index) continue; // P(fn = ft | Gn != Gt) = 0
                       if (ngt != SOMATIC_DIGT::REF)
                       {
                           if (normal_freq_index != ngt)
                               continue;
                           lprior_freq = log_error_mod + ln_cse_rate;
                       }
                       else
                       {
                           if (normal_freq_index == ngt || (consider_norm_contam && normal_freq_index == SOMATIC_DIGT::SIZE))
                               lprior_freq = log_error_mod + ln_one_half;
                           else
                               continue;
                       }
                    }
                    double lsum = lprior_freq + normal_lhood[normal_freq_index] + tumor_lhood[tumor_freq_index];
                    log_sum[index++] = lsum;
                    if(lsum > max_log_sum) max_log_sum = lsum;
                }
            }

            // Calculate log(exp(log_sum[0]-max_log_sum) + ...
            double sum = 0.0;
            for (int i(0); i<index; ++i)
            {
                sum += std::exp(log_sum[i] - max_log_sum);
            }

            // logP(Gn=ngt, Gt=tgt)
            double log_genotype_prior = bare_lnprior[ngt] + ((tgt == 0) ? lnmatch : lnmismatch);

            // log(P(G)) + log(D|G)
            // log(D|G)
            // = log(exp(log_sum[0])+exp(log_sum[1])+...)
            // = max_log_sum + log(exp(log_sum[0]-max_log_sum) + ...
            log_post_prob[ngt][tgt] = log_genotype_prior + max_log_sum + std::log(sum);

            if(log_post_prob[ngt][tgt] > max_log_prob)
            {
                max_log_prob = log_post_prob[ngt][tgt];
                rs.max_gt = DDIGT::get_state(ngt, tgt);
            }
        }
    }

    // Calculate posterior probabilities ( P(G)*P(D|G) )
    double sum_prob = 0.0;
    for (unsigned ngt(0); ngt<SOMATIC_DIGT::SIZE; ++ngt)
    {
        for (unsigned tgt(0); tgt<SOMATIC_STATE::SIZE; ++tgt)
        {
            double prob = std::exp(log_post_prob[ngt][tgt] - max_log_prob); // to prevent underflow
            sum_prob += prob;
        }
    }

    double log_sum_prob = std::log(sum_prob);
    double min_not_somfrom_sum(INFINITY);
    double nonsom_prob = 0.0;

    double post_prob[SOMATIC_DIGT::SIZE][SOMATIC_STATE::SIZE];
    for (unsigned ngt(0); ngt<SOMATIC_DIGT::SIZE; ++ngt)
    {
        double som_prob_given_ngt(0);
        for (unsigned tgt(0); tgt<SOMATIC_STATE::SIZE; ++tgt)
        {
            post_prob[ngt][tgt] = std::exp(log_post_prob[ngt][tgt] - max_log_prob - log_sum_prob);
            if(tgt == 0)    // Non-somatic
            {
                nonsom_prob += post_prob[ngt][tgt];
            }
            else    // Somatic
            {
                som_prob_given_ngt += post_prob[ngt][tgt];
            }
        }

        double err_som_and_ngt = 1.0 - som_prob_given_ngt;
        if (err_som_and_ngt < min_not_somfrom_sum)
        {
            min_not_somfrom_sum=err_som_and_ngt;
            rs.from_ntype_qphred=error_prob_to_qphred(err_som_and_ngt);
            rs.ntype=ngt;
        }
    }
    rs.qphred=error_prob_to_qphred(nonsom_prob);
}
