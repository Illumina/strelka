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

namespace DDIGT
{

void
write_state(const DDIGT::index_t dgt,
            std::ostream& os)
{
    unsigned normal_gt;
    unsigned tumor_gt;
    get_digt_states(dgt,normal_gt,tumor_gt);

    os << SOMATIC_DIGT::label(normal_gt);
    os << "->";
    os << TWO_STATE_SOMATIC::label(tumor_gt);
}

void
write_alt_alleles(unsigned alt_gt,
                  std::ostream& os)
{
    os << id_to_base(alt_gt);
}

}

static const blt_float_t neg_inf = -std::numeric_limits<float>::infinity();

void
calculate_bare_lnprior(const double theta,
        blt_float_t *bare_lnprior)
{
    bare_lnprior[SOMATIC_DIGT::REF]= (blt_float_t) log1p_switch(-(3.*theta)/2.);
    bare_lnprior[SOMATIC_DIGT::HOM]= (blt_float_t) std::log(theta/2.);
    bare_lnprior[SOMATIC_DIGT::HET]= (blt_float_t) std::log(theta);
}

static
float
get_fraction_from_index(int index)
{
    const float ratio_increment(0.5/static_cast<float>(DIGT_GRID::HET_RES+1));
    if (index == SOMATIC_DIGT::REF) return 0.f;
    if (index == SOMATIC_DIGT::HOM) return 1.f;
    if (index == SOMATIC_DIGT::HET) return 0.5f;
    if (index < SOMATIC_DIGT::SIZE+DIGT_GRID::HET_RES) return 1.0f - ratio_increment*(index-SOMATIC_DIGT::SIZE+1);
    return 1.0f - ratio_increment*(index-SOMATIC_DIGT::SIZE+2);
}

inline
int
get_prior_index(
        unsigned ngt,
        unsigned tgt,
        unsigned ft,
        unsigned fn)
{
    return ngt*TWO_STATE_SOMATIC::SIZE*DIGT_GRID::PRESTRAND_SIZE*DIGT_GRID::PRESTRAND_SIZE
            + tgt*DIGT_GRID::PRESTRAND_SIZE*DIGT_GRID::PRESTRAND_SIZE
            + ft*DIGT_GRID::PRESTRAND_SIZE
            + fn;
}

void
set_prior(
        const blt_float_t ssnv_freq_ratio,
        const blt_float_t ln_se_rate,   // ln (shared_error_rate)
        const blt_float_t ln_cse_rate,  // ln (1 - shared_error_rate)
        std::vector<blt_float_t>& ln_freq_given_somatic
        )
{
    blt_float_t log_error_mod = -std::log(static_cast<double>(DIGT_GRID::PRESTRAND_SIZE-1));

    double somatic_prior_normal[DIGT_GRID::PRESTRAND_SIZE] = {};
    somatic_prior_normal[SOMATIC_DIGT::REF] = 0.5;
    somatic_prior_normal[DIGT_GRID::PRESTRAND_SIZE - 1] = 0.5;

    for (unsigned ngt(0); ngt<SOMATIC_DIGT::SIZE; ++ngt)
    {
        // logP(Gn=ngt, Gt=tgt)
        for (unsigned tgt(0); tgt<TWO_STATE_SOMATIC::SIZE; ++tgt) // 0: non-somatic, 1: somatic
        {
            for (unsigned ft(0); ft<DIGT_GRID::PRESTRAND_SIZE; ++ft)
            {
                for (unsigned fn(0); fn<DIGT_GRID::PRESTRAND_SIZE; ++fn)
                {
                    // calculate prior
                    double lprob_f_given_g = 0.0;

                    if(tgt == 0)    // non-somatic
                    {
                       if(fn == ft)
                       {
                           lprob_f_given_g = (fn == ngt) ? ln_cse_rate : ln_se_rate+log_error_mod;
                       }
                       else
                       {
                           lprob_f_given_g = neg_inf;
                       }
                    }
                    else    // somatic
                    {
                       if(fn == ft)
                       {
                           lprob_f_given_g = neg_inf;
                       }
                       else
                       {
                           if (ngt != SOMATIC_DIGT::REF)
                           {
                               lprob_f_given_g = log_error_mod + ((fn == ngt) ? ln_cse_rate : ln_se_rate+log_error_mod);
                           }
                           else
                           {
                               if (get_fraction_from_index(fn) >= ssnv_freq_ratio*get_fraction_from_index(ft))
                                   lprob_f_given_g = neg_inf;
                               else
                               {
                                   lprob_f_given_g = std::log(somatic_prior_normal[fn]) + log_error_mod;
                               }
                           }
                       }
                    }

                    int index = get_prior_index(ngt, tgt, ft, fn);
                    ln_freq_given_somatic[index] = lprob_f_given_g;
                }
            }
        }
    }
}

void
calculate_result_set_grid(
        const blt_float_t* normal_lhood,
        const blt_float_t* tumor_lhood,
        const std::vector<blt_float_t>& ln_somatic_prior,
        const blt_float_t* bare_lnprior,
        const blt_float_t lnmatch,
        const blt_float_t lnmismatch,
        result_set& rs
        )
{
    double log_post_prob[SOMATIC_DIGT::SIZE][TWO_STATE_SOMATIC::SIZE];
    double max_log_prob = neg_inf;

    rs.max_gt = 0;
    for (unsigned ngt(0); ngt<SOMATIC_DIGT::SIZE; ++ngt)
    {
        // logP(Gn=ngt, Gt=tgt)
        double log_diploid_prior_prob = bare_lnprior[ngt];  // logP(Gn)
        for (unsigned tgt(0); tgt<TWO_STATE_SOMATIC::SIZE; ++tgt) // 0: non-somatic, 1: somatic
        {
            double log_prior_prob = log_diploid_prior_prob + ((tgt == 0) ? lnmatch : lnmismatch);
            double max_log_sum = neg_inf;
            double log_sum[DDIGT_GRID::PRESTRAND_SIZE];

            for (unsigned ft(0); ft<DIGT_GRID::PRESTRAND_SIZE; ++ft)
            {
                for (unsigned fn(0); fn<DIGT_GRID::PRESTRAND_SIZE; ++fn)
                {
                    const int prior_index = get_prior_index(ngt, tgt, ft, fn);

                    const unsigned dgt(DDIGT_GRID::get_state(fn, ft));

                    double sum = ln_somatic_prior[prior_index] + normal_lhood[fn] + tumor_lhood[ft];
                    log_sum[dgt] = sum;

                    if(sum > max_log_sum) max_log_sum = sum;
                }
            }

            // calculate log(exp(log_sum[0])+exp(log_sum[1])+...)
            double sum = 0.0;
//            for (int i(0); i<DDIGT_GRID::SIZE; ++i)
            for (int i(0); i<DDIGT_GRID::PRESTRAND_SIZE; ++i)
            {
                sum += std::exp(log_sum[i] - max_log_sum);
            }

            log_post_prob[ngt][tgt] = log_prior_prob + max_log_sum + std::log(sum);

            if(log_post_prob[ngt][tgt] > max_log_prob)
            {
                max_log_prob = log_post_prob[ngt][tgt];
                rs.max_gt = DDIGT::get_state(ngt, tgt);
            }
        }
    }

    // Calculate posterior probabilities
    double sum_prob = 0.0;
    for (unsigned ngt(0); ngt<SOMATIC_DIGT::SIZE; ++ngt)
    {
        for (unsigned tgt(0); tgt<TWO_STATE_SOMATIC::SIZE; ++tgt)
        {
            double prob = std::exp(log_post_prob[ngt][tgt] - max_log_prob); // to prevent underflow
            sum_prob += prob;
        }
    }

    double log_sum_prob = std::log(sum_prob);
    double min_not_somfrom_sum(INFINITY);
    double nonsom_prob = 0.0;

    double post_prob[SOMATIC_DIGT::SIZE][TWO_STATE_SOMATIC::SIZE];
    for (unsigned ngt(0); ngt<SOMATIC_DIGT::SIZE; ++ngt)
    {
        double som_prob_given_ngt(0);
        for (unsigned tgt(0); tgt<TWO_STATE_SOMATIC::SIZE; ++tgt)
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
