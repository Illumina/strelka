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

/// \author Chris Saunders
///

#include "somatic_call_shared.hh"
#include "somatic_indel_grid.hh"

#include "blt_util/math_util.hh"
#include "blt_util/prob_util.hh"
#include "starling_common/indel_util.hh"

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <iostream>
#include <iomanip>
#include <utility>
#include <vector>


#ifdef DEBUG_INDEL_CALL
#include "blt_util/log.hh"
#endif



namespace DDIINDEL_GRID
{

is_nonsom_maker_t::
is_nonsom_maker_t()
    : val(SIZE,false)
{
    for (unsigned gt(0); gt<STAR_DIINDEL_GRID::SIZE; ++gt)
    {
        val[get_state(gt,gt)] = true;
    }
}

const is_nonsom_maker_t is_nonsom;
}



std::ostream&
operator<<(std::ostream& os,
           const DDIINDEL_GRID::index_t dgt)
{

    unsigned normal_gt;
    unsigned tumor_gt;
    DDIINDEL_GRID::get_sdiindel_grid_states(dgt,normal_gt,tumor_gt);

    os << STAR_DIINDEL::label(STAR_DIINDEL_GRID::get_star_diindel_state(normal_gt))
       << "->"
       << STAR_DIINDEL::label(STAR_DIINDEL_GRID::get_star_diindel_state(tumor_gt));

    return os;
}



// For:
//
// homozygous state: S
// frequency grid: G
// number of allele axes: a = 1
// somatic grid size: |G|-a
//
// ln_som_match    = log( 1-P(S) )
// ln_som_mismatch = log( P(S)/(|G|-a) )
//
// :
// indel noise rate: n
// diploid set: D
// non-diploid noise points: N={G\D}
//
// lnp_norm[0..2]   = log( norm[0..2]*(1-n) )
// lnp_norm[3..|G|] = log( n/|N| )
//
//
somatic_indel_caller_grid::
somatic_indel_caller_grid(const strelka_options& opt,
                          const indel_digt_caller& in_caller)
{
    _ln_som_match=(log1p_switch(-opt.somatic_indel_rate));
    _ln_som_mismatch=(std::log(opt.somatic_indel_rate/(static_cast<double>(STAR_DIINDEL_GRID::SIZE-1))));

    std::fill(_bare_lnprior.normal.begin(),_bare_lnprior.normal.end(),0);

    const double* normal_lnprior_genomic(in_caller.lnprior_genomic());

    for (unsigned ngt(0); ngt<STAR_DIINDEL::SIZE; ++ngt)
    {
        _bare_lnprior.normal[ngt] = normal_lnprior_genomic[ngt];
    }
}

static unsigned get_index_from_fraction(float f)
{
    if(f <= 0.5)
    {
        int bin = (int)(f*(STAR_DIINDEL_GRID::HET_RES+1)*2+0.5f);
        if(bin <= 0) return STAR_DIINDEL::NOINDEL;
        else if(bin >= STAR_DIINDEL_GRID::HET_RES+1) return STAR_DIINDEL::HET;
        else return STAR_DIINDEL::SIZE-1+bin;
    }
    else
    {
        int bin = (int)((f-0.5f)*(STAR_DIINDEL_GRID::HET_RES+1)*2+0.5f);
        if(bin <= 0) return STAR_DIINDEL::HET;
        else if(bin >= STAR_DIINDEL_GRID::HET_RES+1) return STAR_DIINDEL::HOM;
        else return STAR_DIINDEL::SIZE-1+STAR_DIINDEL_GRID::HET_RES + bin;
    }
}

static float get_fraction_from_purity(float purity, int gt, int alt_gt)
{
    if(gt == STAR_DIINDEL::NOINDEL)
    {
        if(alt_gt == STAR_DIINDEL::NOINDEL) return 0.f;
        if(alt_gt == STAR_DIINDEL::HET) return (1.0f-purity)/2;
        return 1.0f-purity;
    }
    if(gt == STAR_DIINDEL::HET)
    {
        if(alt_gt == STAR_DIINDEL::NOINDEL) return purity/2;
        if(alt_gt == STAR_DIINDEL::HET) return 0.5f;
        return 1.0f-purity/2;
    }
    // gt == STAR_DIINDEL::HET
    if(alt_gt == STAR_DIINDEL::NOINDEL) return purity;
    if(alt_gt == STAR_DIINDEL::HET) return (1.0f+purity)/2;
    return 1.0f;
}

// prior_normal_lprob = log P(fn | Gn, Gt)
// prior_tumor_lprob = log P(ft | Gn, Gt)
void
somatic_indel_caller_grid::
set_prior(double *prior_normal_lprob,
        double *prior_tumor_lprob,
        const double ref_error_prob
        ) const
{
    const float known_norm_purity = 1.f;
//    const int fn_target = get_index_from_fraction(known_norm_purity);
//    double prob_norm_purity[STAR_DIINDEL_GRID::SIZE];
//    for (unsigned fn(0); fn<STAR_DIINDEL_GRID::SIZE; ++fn)
//    {
//        if(fn == fn_target) prob_norm_purity[fn] = 0.8;
//        if(fn == fn_target+1) prob_norm_purity[fn] = 0.1;
//        if(fn == fn_target-1) prob_norm_purity[fn] = 0.1;
//    }
    const float known_tumor_purity = 0.4f;

    const double prob_miss_normal_fraction = 0.00000001; // 1E-8
    const double prob_miss_tumor_fraction = 0.01;

    double ref_lprob = log1p_switch(-ref_error_prob);  // log(1 - mu)
    double ref_error_lprob = std::log(ref_error_prob);  // log mu
    double noise_lprob = -std::log(static_cast<double>(STAR_DIINDEL_GRID::SIZE-1)); // log P_noise(f)

    double lprob_hit_target_normal = log1p_switch(-prob_miss_normal_fraction); // log (1-prob_miss_normal_fraction)
    double lprob_miss_target_normal = std::log(prob_miss_normal_fraction) - std::log(static_cast<double>(STAR_DIINDEL_GRID::SIZE-1)); // log prob_miss_normal_fraction / (GRID_SIZE-1)

    double lprob_hit_target_tumor = log1p_switch(-prob_miss_tumor_fraction); // log (1-prob_miss_tumor_fraction)
    double lprob_miss_target_tumor = std::log(prob_miss_tumor_fraction) - std::log(static_cast<double>(STAR_DIINDEL_GRID::SIZE-1)); // log prob_miss_tumor_fraction / (GRID_SIZE-1)

    for (unsigned ngt(0); ngt<STAR_DIINDEL::SIZE; ++ngt)
    {
        for (unsigned tgt(0); tgt<STAR_DIINDEL::SIZE; ++tgt)
        {
            float target_norm_frac = get_fraction_from_purity(known_norm_purity, ngt, tgt);
            unsigned target_fn = get_index_from_fraction(target_norm_frac);

            float target_tumor_frac = get_fraction_from_purity(known_tumor_purity, tgt, ngt);
            unsigned target_tn = get_index_from_fraction(target_tumor_frac);

            int base_index = (STAR_DIINDEL::SIZE*ngt + tgt)*STAR_DIINDEL_GRID::SIZE;
            for(unsigned f(0); f<STAR_DIINDEL_GRID::SIZE; f++)
            {
                int index = base_index + f;
                if(ngt == tgt)  // non-somatic
                {
                    if (f==ngt)
                    {
                        prior_normal_lprob[index] = ref_lprob;
                        prior_tumor_lprob[index] = ref_lprob;
                    }
                    else
                    {
                        prior_normal_lprob[index] = ref_error_lprob + noise_lprob;
                        prior_tumor_lprob[index] = ref_error_lprob + noise_lprob;
                    }
                }
                else    // somatic state
                {
                    if (f==target_fn) prior_normal_lprob[index] = lprob_hit_target_normal;
                    else prior_normal_lprob[index] = lprob_miss_target_normal;

                    if (f==target_tn) prior_tumor_lprob[index] = lprob_hit_target_tumor;
                    else prior_tumor_lprob[index] = lprob_miss_target_tumor;
                }
            }
        }
    }
}


typedef somatic_indel_call::result_set result_set;


static
void
get_indel_het_grid_lhood(const starling_base_options& opt,
                         const starling_base_deriv_options& dopt,
                         const starling_sample_options& sample_opt,
                         const double indel_error_lnp,
                         const double indel_real_lnp,
                         const double ref_error_lnp,
                         const double ref_real_lnp,
                         const indel_key& ik,
                         const indel_data& id,
                         const bool is_include_tier2,
                         const bool is_use_alt_indel,
                         double* const lhood)
{
    static const unsigned lsize(STAR_DIINDEL_GRID::HET_RES*2);
    for (unsigned gt(0); gt<(lsize); ++gt) lhood[gt] = 0.;

    static const double ratio_increment(0.5/static_cast<double>(STAR_DIINDEL_GRID::HET_RES+1));
    for (unsigned i(0); i<STAR_DIINDEL_GRID::HET_RES; ++i)
    {
        const double het_ratio((i+1)*ratio_increment);
        indel_digt_caller::get_high_low_het_ratio_lhood(opt,dopt,
                                                        sample_opt,
                                                        indel_error_lnp,indel_real_lnp,
                                                        ref_error_lnp,ref_real_lnp,
                                                        ik,id,het_ratio,
                                                        is_include_tier2,is_use_alt_indel,
                                                        lhood[i],
                                                        lhood[lsize-(i+1)]);
    }
}

static
double
get_lprob_f_given_ngt_tgt(
        const double *prior_lprob,
        const unsigned f,
        const unsigned ngt,
        const unsigned tgt)
{
    int index = (STAR_DIINDEL::SIZE*ngt + tgt)*STAR_DIINDEL_GRID::SIZE + f;
    return prior_lprob[index];
}

void
somatic_indel_caller_grid::calculate_result_set(
    const double *prior_normal_lprob,
    const double *prior_tumor_lprob,
    const double ref_error_lnp,
    const double ref_real_lnp,
    const double* normal_lhood,
    const double* tumor_lhood,
    result_set& rs) const
{
    double log_post_prob[STAR_DIINDEL::SIZE][STAR_DIINDEL::SIZE];
    double max_log_prob = -INFINITY;
    rs.max_gt=0;

    double log_one_third = std::log(1.0/3.0);
    double log_two_third = std::log(2.0/3.0);
    double log_half = std::log(0.5);

    for (int ngt(0); ngt<STAR_DIINDEL::SIZE; ++ngt)
    {
        // logP(Gn=ngt, Gt=tgt)
        double log_prior_prob = _bare_lnprior.normal[ngt];

        for (int tgt(0); tgt<STAR_DIINDEL::SIZE; ++tgt)
        {
            if(ngt == tgt)
            {
                log_prior_prob += ref_real_lnp;     // log(1-P(S))
            }
            else
            {
                log_prior_prob += ref_error_lnp;
                if(ngt == STAR_DIINDEL::NOINDEL)
                {
                    if(tgt == STAR_DIINDEL::HET) log_prior_prob += log_two_third;   // NOINDEL, HET
                    else log_prior_prob += log_one_third;   // NOINDEL, HOM
                }
                else if(ngt == STAR_DIINDEL::HET)
                {
                    log_prior_prob += log_half;   // HET, NOINDEL or HET, HOM
                }
                else    // ngt == STAR_DIINDEL::HOM
                {
                    if(tgt == STAR_DIINDEL::NOINDEL) log_prior_prob += log_one_third;   // HOM, NOINDEL
                    else log_prior_prob += log_two_third;   // HOM, HET
                }
            }

            log_post_prob[ngt][tgt] = log_prior_prob;

            double max_log_sum = -INFINITY;
            double log_sum[STAR_DIINDEL_GRID::SIZE*STAR_DIINDEL_GRID::SIZE];
            for (unsigned fn(0); fn<STAR_DIINDEL_GRID::SIZE; ++fn)
            {
                for (unsigned ft(0); ft<STAR_DIINDEL_GRID::SIZE; ++ft)
                {
                    int i = fn*STAR_DIINDEL_GRID::SIZE + ft;
                    // logP(fn|Gn=ngt,Gt=tgt)
                    double lprob_fn_given_gn_gt = get_lprob_f_given_ngt_tgt(prior_normal_lprob, fn, ngt, tgt);
                    // logP(ft|Gn=ngt,Gt=tgt)
                    double lprob_ft_given_gn_gt = get_lprob_f_given_ngt_tgt(prior_tumor_lprob, ft, ngt, tgt);

                    // logP(fn|Gn=ngt,Gt=tgt)+logP(ft|Gn=ngt,Gt=tgt)+logP(Dn|fn)+logP(Dt|ft)
                    double sum = lprob_fn_given_gn_gt + lprob_ft_given_gn_gt + normal_lhood[fn] + tumor_lhood[ft];
                    log_sum[i] = sum;

                    if(sum > max_log_sum) max_log_sum = sum;
                }
            }

            // calculate log(exp(log_sum[0])+exp(log_sum[1])+...)
            double sum = 0.0;
            for (int i(0); i<STAR_DIINDEL_GRID::SIZE*STAR_DIINDEL_GRID::SIZE; ++i)
            {
                sum += std::exp(log_sum[i] - max_log_sum);
            }

            log_post_prob[ngt][tgt] += max_log_sum + std::log(sum);

            if(log_post_prob[ngt][tgt] > max_log_prob)
            {
                max_log_prob = log_post_prob[ngt][tgt];
                rs.max_gt = DDIINDEL_GRID::get_state(ngt,tgt);
            }
        }
    }

    // log prob to prob, calculate sum
    double post_prob[STAR_DIINDEL::SIZE][STAR_DIINDEL::SIZE];
    double sum_prob = 0.0;
    for (unsigned ngt(0); ngt<STAR_DIINDEL::SIZE; ++ngt)
    {
        for (unsigned tgt(0); tgt<STAR_DIINDEL::SIZE; ++tgt)
        {
            double prob = std::exp(log_post_prob[ngt][tgt] - max_log_prob); // to prevent underflow
            sum_prob += prob;
//            post_prob[ngt][tgt] = prob;
        }
    }

//    rs.sindel_from_ntype_qphred=error_prob_to_qphred(not_somfrom_sum);
//    rs.ntype=sgt;

    // normalize probabilities, determine QSI, NT and QSI_NT
    double log_sum_prob = std::log(sum_prob);
    double min_not_somfrom_sum(0);
    double nonsomatic_sum(0);
    for (unsigned ngt(0); ngt<STAR_DIINDEL::SIZE; ++ngt)
    {
        double som_prob = 0.0;
        for (unsigned tgt(0); tgt<STAR_DIINDEL::SIZE; ++tgt)
        {
//            post_prob[ngt][tgt] /= sum_prob;
            post_prob[ngt][tgt] = std::exp(log_post_prob[ngt][tgt] - max_log_prob - log_sum_prob);
            if(tgt != ngt)
            {
                som_prob += post_prob[ngt][tgt];
            }
            else
            {
                nonsomatic_sum += post_prob[ngt][tgt];
            }
        }

        double not_somfrom_sum = 1.0 - som_prob;
        if ((ngt==0) || (not_somfrom_sum<min_not_somfrom_sum))
        {
            min_not_somfrom_sum=not_somfrom_sum;
            rs.sindel_from_ntype_qphred=error_prob_to_qphred(not_somfrom_sum);
            rs.ntype=ngt;
        }
    }
    rs.sindel_qphred=error_prob_to_qphred(nonsomatic_sum);

    if (0 == rs.sindel_qphred) return;
}


#ifdef SOMATIC_DEBUG
static
void
debug_dump_indel_lhood(const double* lhood,
                       std::ostream& os)
{

    double pprob[STAR_DIINDEL_GRID::SIZE];
    for (unsigned gt(0); gt<STAR_DIINDEL_GRID::SIZE; ++gt)
    {
        pprob[gt] = lhood[gt];
    }

    unsigned max_gt(0);
    normalize_ln_distro(pprob,pprob+STAR_DIINDEL_GRID::SIZE,max_gt);

    os << std::setprecision(3) << std::fixed;
    for (unsigned gt(0); gt<STAR_DIINDEL_GRID::SIZE; ++gt)
    {
        os << static_cast<STAR_DIINDEL::index_t>(STAR_DIINDEL_GRID::get_star_diindel_state(gt)) << ": " << -std::log(pprob[gt]) << " ";
        os << "\n";
    }
    os.unsetf(std::ios::fixed);
}
#endif



static
bool
is_multi_indel_allele(const starling_base_deriv_options& dopt,
                      const indel_data& normal_id,
                      const indel_data& tumor_id,
                      const bool is_include_tier2,
                      bool& is_overlap)
{
    static const bool is_use_alt_indel(true);
    static const double min_explained_count_fraction(.9);

    enum iallele_t
    {
        INDEL = -2,
        REF = -1
    };

    // get total pprob:
    read_path_scores total_pprob;
    get_sum_path_pprob(dopt,normal_id,is_include_tier2,is_use_alt_indel,total_pprob,true);
    get_sum_path_pprob(dopt,tumor_id,is_include_tier2,is_use_alt_indel,total_pprob,false);

    // next determine the top two indel alleles:
    std::vector<std::pair<double,int> > scores;
    scores.push_back(std::make_pair(-total_pprob.indel,static_cast<int>(INDEL)));
    scores.push_back(std::make_pair(-total_pprob.ref,static_cast<int>(REF)));
    const read_path_scores::alt_indel_t& ai(total_pprob.alt_indel);
    const int ais(ai.size());
    for (int i(0); i<ais; ++i)
    {
        scores.push_back(std::make_pair(-ai[i].second,i));
    }

    sort(scores.begin(),scores.end());

#if 0
    std::cerr << "BUG: normal_id: " << normal_id;
    std::cerr << "BUG: tumor_id: " << tumor_id;
    std::cerr << "BUG: total_pprob: " << total_pprob << "\n";
    std::cerr << "BUG: max1_id,max2_id " << scores[0].second << " " << scores[1].second << "\n";
#endif

    // If the top two alleles are both alternate indels, check that
    // they interfere with each other.  If not, we are forced to make
    // the conservative assumption that they occur as part of the same
    // haplotype:
    //
    assert(scores.size() >= 2);
    while (scores[0].second>=0 && scores[1].second>=0)
    {
        if (is_indel_conflict(ai[scores[0].second].first,ai[scores[1].second].first))
        {
            break;
        }
        scores.erase(scores.begin()+1);
        assert(scores.size() >= 2);
    }

    if ((scores[0].second!=INDEL) && (scores[1].second!=INDEL)) return true;
    if (scores.size() >= 3)
    {
        const double top_prob(scores[0].first+scores[1].first);
        const double top_frac(top_prob/(top_prob+scores[2].first));
        if (top_frac<min_explained_count_fraction) return true;
    }


    // the rejection criteria is resolved at this point, but one more
    // flag is set below as an interesting utility for users to
    // quickly find/filter the 'legitimate' overlapping indels:
    //
    is_overlap=((scores[0].second!=REF) && (scores[1].second!=REF));

    return false;
}



///
///
void
somatic_indel_caller_grid::
get_somatic_indel(const strelka_options& opt,
                  const strelka_deriv_options& dopt,
                  const starling_sample_options& normal_opt,
                  const starling_sample_options& tumor_opt,
                  const double indel_error_prob,
                  const double ref_error_prob,
                  const indel_key& ik,
                  const indel_data& normal_id,
                  const indel_data& tumor_id,
                  const bool is_use_alt_indel,
                  somatic_indel_call& sindel) const
{
    // for now, lhood calculation of tumor and normal are independent:

    // get likelihood of each genotype
    static const bool is_normal_het_bias(false);
    static const bool is_tumor_het_bias(false);
    static const double normal_het_bias(0.0);
    static const double tumor_het_bias(0.0);
    double normal_lhood[STAR_DIINDEL_GRID::SIZE];
    double tumor_lhood[STAR_DIINDEL_GRID::SIZE];

    const int array_size = STAR_DIINDEL::SIZE*STAR_DIINDEL::SIZE*STAR_DIINDEL_GRID::SIZE;
    double prior_normal_lprob[array_size];
    double prior_tumor_lprob[array_size];

    sindel.is_forced_output=(normal_id.is_forced_output || tumor_id.is_forced_output);

    const double indel_error_lnp(std::log(indel_error_prob));
    const double indel_real_lnp(std::log(1.-indel_error_prob));
    const double ref_error_lnp(std::log(ref_error_prob));
    const double ref_real_lnp(std::log(1.-ref_error_prob));

    static const unsigned n_tier(2);
    std::array<result_set,n_tier> tier_rs;
    for (unsigned i(0); i<n_tier; ++i)
    {
        const bool is_include_tier2(i==1);
        if (is_include_tier2)
        {
            if (! opt.tier2.is_tier2()) continue;
            if (tier_rs[0].sindel_qphred==0)
            {
                if (! sindel.is_forced_output)   // if forced output then there's still a point to computing tier2
                {
                    tier_rs[1].sindel_qphred=0;
                    continue;
                }
            }
        }

        static const bool is_somatic_multi_indel_filter(true);
#if 0
        std::cerr << "BUG: testing tier/ik: " << i << " " << ik;
#endif
        if (is_somatic_multi_indel_filter)
        {
            const bool ismulti(is_multi_indel_allele(dopt,normal_id,tumor_id,is_include_tier2,tier_rs[i].is_overlap));
            if (ismulti)
            {
                tier_rs[i].sindel_qphred=0;
#if 0
                std::cerr << "BUG: rejected\n";
#endif
                continue;
            }
        }

        indel_digt_caller::get_indel_digt_lhood(opt,dopt,normal_opt,
                                                indel_error_prob,ref_error_prob,ik,normal_id,
                                                is_normal_het_bias,normal_het_bias,
                                                is_include_tier2,is_use_alt_indel,
                                                normal_lhood);
        indel_digt_caller::get_indel_digt_lhood(opt,dopt,tumor_opt,
                                                indel_error_prob,ref_error_prob,ik,tumor_id,
                                                is_tumor_het_bias,tumor_het_bias,
                                                is_include_tier2,is_use_alt_indel,
                                                tumor_lhood);

        get_indel_het_grid_lhood(opt,dopt,normal_opt,
                                 indel_error_lnp,indel_real_lnp,
                                 ref_error_lnp,ref_real_lnp,
                                 ik,normal_id,
                                 is_include_tier2,is_use_alt_indel,
                                 normal_lhood+STAR_DIINDEL::SIZE);
        get_indel_het_grid_lhood(opt,dopt,tumor_opt,
                                 indel_error_lnp,indel_real_lnp,
                                 ref_error_lnp,ref_real_lnp,
                                 ik,tumor_id,
                                 is_include_tier2,is_use_alt_indel,
                                 tumor_lhood+STAR_DIINDEL::SIZE);

        set_prior(prior_normal_lprob, prior_tumor_lprob,
                ref_error_prob);

        calculate_result_set(prior_normal_lprob, prior_tumor_lprob,
                ref_error_lnp,ref_real_lnp,
                normal_lhood,tumor_lhood,tier_rs[i]);
    }

    if (! sindel.is_forced_output)
    {
        if (tier_rs[0].sindel_qphred==0 ||
            tier_rs[1].sindel_qphred==0) return;
    }

    sindel.sindel_tier=0;
    if (opt.tier2.is_tier2())
    {
        if (tier_rs[0].sindel_qphred > tier_rs[1].sindel_qphred)
        {
            sindel.sindel_tier=1;
        }
    }

    sindel.sindel_from_ntype_tier=0;
    if (opt.tier2.is_tier2())
    {
        if (tier_rs[0].sindel_from_ntype_qphred > tier_rs[1].sindel_from_ntype_qphred)
        {
            sindel.sindel_from_ntype_tier=1;
        }
    }

    sindel.rs=tier_rs[sindel.sindel_from_ntype_tier];

    if (tier_rs[0].ntype != tier_rs[1].ntype)
    {
        // catch NTYPE conflict states:
        sindel.rs.ntype = NTYPE::CONFLICT;
        sindel.rs.sindel_from_ntype_qphred = 0;

    }
    else
    {
        // classify NTYPE:
        //

        // convert diploid genotype into more limited ntype set:
        //
        if       (sindel.rs.ntype==STAR_DIINDEL::NOINDEL)
        {
            sindel.rs.ntype=NTYPE::REF;
        }
        else if (sindel.rs.ntype==STAR_DIINDEL::HOM)
        {
            sindel.rs.ntype=NTYPE::HOM;
        }
        else
        {
            sindel.rs.ntype=NTYPE::HET;
#if 0
        }
        else if (sindel.rs.ntype==STAR_DIINDEL::HET)
        {
            sindel.rs.ntype=NTYPE::HET;
        }
        else
        {
            assert(0);
#endif
        }
    }

    sindel.rs.sindel_qphred = tier_rs[sindel.sindel_tier].sindel_qphred;
}
