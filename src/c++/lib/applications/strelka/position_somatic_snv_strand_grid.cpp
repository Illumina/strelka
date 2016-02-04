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

#include "position_somatic_snv_strand_grid.hh"
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

//#define SOMATIC_DEBUG

#ifdef SOMATIC_DEBUG
#include <iostream>
#include <iomanip>
#endif


static const blt_float_t one_third(1./3.);
static const blt_float_t ln_one_third(std::log(one_third));
static const blt_float_t one_half(1./2.);
static const blt_float_t ln_one_half(std::log(one_half));



// compute just the non-strand-bias portion of the normal marginal
// prior given p(signal), p(no-strand noise), p(strand-bias noise)
//
static
void
get_nostrand_marginal_prior(const blt_float_t* normal_lnprior,
                            const blt_float_t sse_rate,
                            const blt_float_t sseb_fraction,
                            std::vector<blt_float_t>& grid_normal_lnprior)
{
    const blt_float_t strand_sse_rate(sse_rate*sseb_fraction);
    const blt_float_t nostrand_sse_rate(sse_rate-strand_sse_rate);

    const blt_float_t ln_csse_rate( log1p_switch(-sse_rate) );
    //    const blt_float_t ln_strand_sse_rate( std::log(strand_sse_rate) );
    const blt_float_t ln_nostrand_sse_rate( std::log(nostrand_sse_rate) );

    // fill in normal sample prior for canonical diploid allele frequencies:
    for (unsigned ngt(0); ngt<DIGT_SIMPLE::SIZE; ++ngt)
    {
        grid_normal_lnprior[ngt] = (normal_lnprior[ngt]+ln_csse_rate);
    }

    static const blt_float_t error_mod( -std::log(static_cast<blt_float_t>(DIGT_SGRID::HET_RES*2)) );
    // fill in normal sample prior for 'noise' frequencies:
    for (unsigned ngt(DIGT_SIMPLE::SIZE); ngt<DIGT_SGRID::PRESTRAND_SIZE; ++ngt)
    {
        grid_normal_lnprior[ngt] = (normal_lnprior[DIGT_SIMPLE::HET]+ln_nostrand_sse_rate+error_mod);
    }
}



//
// The method only pre-computes the marginal normal allele-frequency
// component of the prior, the prior is expanded to include the tumor
// allele frequency during the posterior computation
//
//
// For:
// somatic state S
// snv noise rate n
// snv strand-biased noise fraction x
//
// ln_csse_rate = log( 1-n )
// ln_strand_sse_rate = log( nx )
// ln_nostrand_sse_rate = log( n(1-x) )
//
// grid_normal_lnprior is the mixture of normal diploid probabilities with uniform noise
//
static
void
get_prior(const blt_float_t* normal_lnprior,
          const blt_float_t sse_rate,
          const blt_float_t sseb_fraction,
          const blt_float_t somatic_normal_noise_rate,
          const bool is_somatic_normal_noise_rate,
          std::vector<blt_float_t>& grid_normal_lnprior,
          std::vector<blt_float_t>& somatic_marginal_lnprior)
{
    get_nostrand_marginal_prior(normal_lnprior,sse_rate,sseb_fraction,grid_normal_lnprior);
    if (is_somatic_normal_noise_rate)
    {
        get_nostrand_marginal_prior(normal_lnprior,somatic_normal_noise_rate,0,somatic_marginal_lnprior);
    }
    else
    {
        get_nostrand_marginal_prior(normal_lnprior,sse_rate,sseb_fraction,somatic_marginal_lnprior);
    }

    const blt_float_t strand_sse_rate(sse_rate*sseb_fraction);
    //    const blt_float_t nostrand_sse_rate(sse_rate-strand_sse_rate);

    //    const blt_float_t ln_csse_rate( log1p_switch(-sse_rate) );
    const blt_float_t ln_strand_sse_rate( std::log(strand_sse_rate) );
    //    const blt_float_t ln_nostrand_sse_rate( std::log(nostrand_sse_rate) );

    static const blt_float_t error_mod( -std::log(static_cast<blt_float_t>(DIGT_SGRID::HET_RES*2)) );

    // flaw in using error mod -- because strand state could exist and
    // be detectable at the canonical gt frequencies, we leave it the
    // same with the goal of stable performance as the user changes
    // the fraction term on the command-line:
    //
    // TODO: unclear comment?
    //

    for (unsigned ngt(DIGT_SGRID::PRESTRAND_SIZE); ngt<(DIGT_SGRID::SIZE); ++ngt)
    {
        grid_normal_lnprior[ngt] = (grid_normal_lnprior[DIGT_SIMPLE::HET]+ln_strand_sse_rate+error_mod);
    }

#ifdef SOMATIC_DEBUG
    double throwaway_sum(0);
    for (unsigned sgt(0); sgt<n_strand_het_axes; ++sgt)
    {
        throwaway_sum += std::exp(strand_axis_prior[sgt]);
    }

    throwaway_sum /= 2.; // we're only using half-axes for the states represented
    throwaway_sum = 1. - throwaway_sum;

    throwaway_sum *= std::exp(ln_strand_sse_rate);

    check_ln_distro(grid_normal_lnprior.begin(),
                    grid_normal_lnprior.end(),
                    "somatic prior",
                    0.0001,
                    1.-throwaway_sum);
#endif
}



somatic_snv_caller_strand_grid::
somatic_snv_caller_strand_grid(const strelka_options& opt,
                               const pprob_digt_caller& pd_caller)
    : _opt(opt)
{
    _ln_som_match=(log1p_switch(-opt.somatic_snv_rate));
//    _ln_som_mismatch=(std::log(opt.somatic_snv_rate/(static_cast<blt_float_t>((DIGT_SGRID::PRESTRAND_SIZE)-1))));
    _ln_som_mismatch=(std::log(opt.somatic_snv_rate));

    const blt_float_t strand_sse_rate(opt.shared_site_error_rate*opt.shared_site_error_strand_bias_fraction);
    const blt_float_t nostrand_sse_rate(opt.shared_site_error_rate-strand_sse_rate);

    _ln_csse_rate = log1p_switch(-_opt.shared_site_error_rate);
    _ln_nostrand_sse_rate = std::log(nostrand_sse_rate);
    _bare_lnprior = pd_caller.lnprior_genomic();

    prior_set& ps(_lnprior);
    std::fill(ps.normal.begin(),ps.normal.end(),0);
    std::fill(ps.normal_poly.begin(),ps.normal_poly.end(),0);

    get_prior(pd_caller.lnprior_genomic(),
              opt.shared_site_error_rate,
              opt.shared_site_error_strand_bias_fraction,
              opt.site_somatic_normal_noise_rate,
              opt.is_site_somatic_normal_noise_rate,
              ps.normal,
              ps.somatic_marginal);
    get_prior(pd_caller.lnprior_polymorphic(),
              opt.shared_site_error_rate,
              opt.shared_site_error_strand_bias_fraction,
              opt.site_somatic_normal_noise_rate,
              opt.is_site_somatic_normal_noise_rate,
              ps.normal_poly,
              ps.somatic_marginal_poly);

    // special nostrand distro is used for somatic_gvcf:
    get_nostrand_marginal_prior(pd_caller.lnprior_genomic(),opt.shared_site_error_rate,0,ps.normal_nostrand);
    get_nostrand_marginal_prior(pd_caller.lnprior_polymorphic(),opt.shared_site_error_rate,0,ps.normal_poly_nostrand);
}



// calculate probability of strand-specific noise
//
// accelerated version with no hyrax q-val mods:
//
// the ratio key can be used as a proxy for the het ratio to look up cached results:
//
static
void
get_strand_ratio_lhood_spi(
    const snp_pos_info& pi,
    const unsigned ref_gt,
    const blt_float_t het_ratio,
    const unsigned het_ratio_index,
    het_ratio_cache<2>& hrcache,
    blt_float_t* lhood)
{
    // het_ratio is the expected allele frequency of noise on the
    // noise-strand, or "on-strand" below. All possible ratio values
    // (there should be very few), have an associated index value for
    // caching.
    //
    const blt_float_t chet_ratio(1.-het_ratio);

    // In this situation every basecall falls into 1 of 4 states:
    //
    // 0: off-strand non-reference allele (0)
    // 1: on-strand non-reference allele (het_ratio)  (cached)
    // 2: on-strand agrees with the reference (chet_ratio) (cached)
    // 3: off-strand agree with the reference (1)
    //
    // The off-strand states are not cached below because they're
    // simpler to compute
    //
    blt_float_t lhood_fwd[DIGT_SGRID::STRAND_SIZE]; // "on-strand" is fwd
    blt_float_t lhood_rev[DIGT_SGRID::STRAND_SIZE]; // "on-strand" is rev

    for (unsigned i(0); i<DIGT_SGRID::STRAND_SIZE; ++i)
    {
        lhood_fwd[i] = 0;
        lhood_rev[i] = 0;
    }

    static const unsigned n_strand_het_axes(DIGT_SGRID::STRAND_SIZE);

    for (const base_call& bc : pi.calls)
    {
        std::pair<bool,cache_val<2>*> ret(hrcache.get_val(bc.get_qscore(),het_ratio_index));
        cache_val<2>& cv(*ret.second);

        // compute results only if they aren't already cached:
        //
        if (! ret.first)
        {
            const blt_float_t eprob(bc.error_prob());
            const blt_float_t ceprob(1.-eprob);
            // cached value [0] refers to state 2 above: on-strand
            // reference allele
            cv.val[0]=(std::log((ceprob)*chet_ratio+((eprob)*one_third)*het_ratio));
            // cached value [1] refers to state 1 above: on-strand
            // non-reference allele
            cv.val[1]=(std::log((ceprob)*het_ratio+((eprob)*one_third)*chet_ratio));
        }

        const uint8_t obs_id(bc.base_id);

        if (obs_id==ref_gt)
        {
            const blt_float_t val_off_strand(bc.ln_comp_error_prob());
            const blt_float_t val_fwd(bc.is_fwd_strand ? cv.val[0] : val_off_strand);
            const blt_float_t val_rev(bc.is_fwd_strand ? val_off_strand : cv.val[0]);
            for (unsigned sgt(0); sgt<n_strand_het_axes; ++sgt)
            {
                lhood_fwd[sgt] += val_fwd;
                lhood_rev[sgt] += val_rev;
            }
        }
        else
        {
            const blt_float_t val_off_strand(bc.ln_error_prob()+ln_one_third);
            const blt_float_t val_fwd(bc.is_fwd_strand ? cv.val[1] : val_off_strand);
            const blt_float_t val_rev(bc.is_fwd_strand ? val_off_strand : cv.val[1]);

            const unsigned match_strand_state(obs_id>ref_gt ? obs_id-1 : obs_id);
            for (unsigned sgt(0); sgt<n_strand_het_axes; ++sgt)
            {
                if (sgt==match_strand_state)
                {
                    lhood_fwd[sgt] += val_fwd;
                    lhood_rev[sgt] += val_rev;
                }
                else
                {
                    lhood_fwd[sgt] += val_off_strand;
                    lhood_rev[sgt] += val_off_strand;
                }
            }
        }
    }

    for (unsigned i(0); i<DIGT_SGRID::STRAND_SIZE; ++i)
    {
        lhood[i] = log_sum(lhood_fwd[i],lhood_rev[i])+ln_one_half;
    }
}



// Fill in the noise portions of the likelihood function for the
// regions where we expect strand bias noise (a minor allele frequency
// < 0.5 + reference allele):
//
static
void
get_diploid_strand_grid_lhood_spi(
    const snp_pos_info& pi,
    const unsigned ref_gt,
    blt_float_t* const lhood)
{
    // TODO: make this thread safe...
    //
    static het_ratio_cache<2> hrcache;

    // get likelihood of each genotype
    for (unsigned gt(0); gt<(DIGT_SGRID::STRAND_STATE_SIZE); ++gt) lhood[gt] = 0.;

    static const blt_float_t ratio_increment(0.5/static_cast<blt_float_t>(DIGT_SGRID::HET_RES+1));
    for (unsigned i(0); i<DIGT_SGRID::HET_RES; ++i)
    {
        const blt_float_t het_ratio((i+1)*ratio_increment);
        get_strand_ratio_lhood_spi(pi,ref_gt,het_ratio,i,hrcache,
                                   lhood+i);
    }
}



typedef somatic_snv_genotype_grid::result_set result_set;



#ifdef SOMATIC_DEBUG
#if 0
static
void
debug_dump_ddigt_lhood(const blt_float_t* lhood,
                       std::ostream& os)
{
    double pprob[DDIGT::SIZE]; //intentionally run at higher float-resolution
    for (unsigned gt(0); gt<DDIGT::SIZE; ++gt)
    {
        pprob[gt] = lhood[gt];
    }

    unsigned max_gt(0);
    normalize_ln_distro(pprob,pprob+DDIGT::SIZE,max_gt);

    os << std::setprecision(3) << std::fixed;
    for (unsigned ngt(0); ngt<DIGT::SIZE; ++ngt)
    {
        for (unsigned tgt(0); tgt<DIGT::SIZE; ++tgt)
        {
            const unsigned dgt(DDIGT::get_state(ngt,tgt));
            os << static_cast<DDIGT::index_t>(dgt) << ": " << -std::log(pprob[dgt]) << " ";
        }
        os << "\n";
    }
    os.unsetf(std::ios::fixed);
}
#endif



static
void
sort_n_dump(const std::string& label,
            std::vector<double>& distro,
            std::vector<double>& distro2,
            const unsigned ref_gt)
{
    static const unsigned topn(25);
    std::ostream& os(log_os);

    unsigned max_gt(0);
    normalize_ln_distro(distro.begin(),distro.end(),max_gt);
    unsigned max_gt2(0);
    normalize_ln_distro(distro2.begin(),distro2.end(),max_gt2);

    const unsigned ds(distro.size());
    assert(ds == distro2.size());
    std::vector<std::pair<std::pair<double,double>,unsigned> > idistro(ds);
    for (unsigned i(0); i<ds; ++i)
    {
        idistro[i] = std::make_pair(std::make_pair(-distro[i],-distro2[i]),i);
    }

    const unsigned topn_used(std::min(topn,ds));
    std::partial_sort(idistro.begin(),
                      idistro.begin()+topn_used,
                      idistro.end());

    os << std::setprecision(3) << std::fixed;
    for (unsigned i(0); i<topn_used; ++i)
    {
        os << label << ": " << i << " ";
        DDIGT_SGRID::write_full_state(static_cast<DDIGT_SGRID::index_t>(idistro[i].second),
                                      ref_gt,os);
        os << " " << -std::log(-idistro[i].first.first)
           << " " << -std::log(-idistro[i].first.second) << "\n";
    }
    os.unsetf(std::ios::fixed);
}


#endif



#if 0
static
float
nonsomatic_gvcf_prior(
    const somatic_snv_caller_strand_grid::prior_set& /*pset*/,
    const unsigned ngt,
    const unsigned tgt)
{
    static const float ln_ref_som_match=std::log(0.5);
    static const float ln_ref_som_mismatch=(std::log(0.5/static_cast<blt_float_t>(DIGT_SGRID::PRESTRAND_SIZE-1)));

    if (tgt==ngt)
    {
        return /*pset.normal_poly_nostrand[ngt]+*/ln_ref_som_match;
    }
    else
    {
        return /*pset.normal_poly_nostrand[ngt]+*/ln_ref_som_mismatch;
    }
}
#endif



/// expanded nonsomatic definition for the purpose of somatic gvcf output:
#if 0
static
bool
is_gvcf_nonsomatic_state(
    const unsigned ngt,
    const unsigned tgt)
{
    return (ngt==tgt);
}
#endif


/// expanded definition of 'nonsomatic' for the purpose of providing somatic gVCF output:
//static
//float
//gvcf_nonsomatic_gvcf_prior(
//    const unsigned ngt,
//    const unsigned tgt)
//{
//    if (ngt == tgt)
//    {
//        static const float lone(std::log(1.f));
//        return lone;
//    }
//    else if (ngt<N_BASE)
//    {
//        static const float lhalf(std::log(0.5f));
//        return lhalf;
//    }
//    else
//    {
//        static const float lzero(-std::numeric_limits<float>::infinity());
//        return lzero;
//    }
//}

float get_fraction_from_index_snv(int index)
{
    const float ratio_increment(0.5/static_cast<float>(DIGT_SGRID::HET_RES+1));
    if(index == DIGT_SIMPLE::REF) return 0.f;
    if(index == DIGT_SIMPLE::HOM) return 1.f;
    if(index == DIGT_SIMPLE::HET) return 0.5f;
    if(index < DIGT_SIMPLE::SIZE+DIGT_SGRID::HET_RES) return 1.0f - ratio_increment*(index-DIGT_SIMPLE::SIZE+1);
    return 1.0f - ratio_increment*(index-DIGT_SIMPLE::SIZE+2);
}

static
void
calculate_result_set_grid(
        const bool isComputeNonSomatic,
        const blt_float_t* normal_lhood,
        const blt_float_t* tumor_lhood,
        const blt_float_t* bare_lnprior_normal,
        const blt_float_t ln_se_rate,   // ln (somatic error rate)
        const blt_float_t ln_cse_rate,  // ln (1 - somatic_error_rate)
        const blt_float_t lnmatch,
        const blt_float_t lnmismatch,
        const bool is_forced_output,
        result_set& rs
        )
{
    // Debugging
//    for(int i=0; i<DIGT_SGRID::PRESTRAND_SIZE; ++i)
//    {
//        printf("%f\t%lf\t%lf\n", get_fraction_from_index_snv(i), normal_lhood[i], tumor_lhood[i]);
//    }

    bool is_normal_contaminated = true;
    // Calculate posterior probabilities
    double log_post_prob[DIGT_SIMPLE::SIZE][TWO_STATE_SOMATIC::SIZE];
    double max_log_prob = -INFINITY;

    blt_float_t log_error_mod = -std::log(static_cast<double>(DIGT_SGRID::PRESTRAND_SIZE-1));

    // TODO: precompute this
    double somatic_prior_normal[DIGT_SGRID::PRESTRAND_SIZE] = {};
    somatic_prior_normal[DIGT_SIMPLE::REF] = 0.5; // fn = 0.0
    somatic_prior_normal[DIGT_SGRID::PRESTRAND_SIZE - 1] = 0.5;  // fn = 0.05
//    somatic_prior_normal[DIGT_SGRID::PRESTRAND_SIZE - 2] = 0.05;  // fn = 0.1

//    double somatic_prior[DIGT_SGRID::PRESTRAND_SIZE][DIGT_SGRID::PRESTRAND_SIZE] = {};
//    for (unsigned ft(0); ft<DIGT_SGRID::PRESTRAND_SIZE; ++ft)
//    {
//        for (unsigned fn(0); fn<DIGT_SGRID::PRESTRAND_SIZE; ++fn)
//        {
//            somatic_prior[std::log(somatic_prior_normal[fn]) + std::log(somatic_prior_tumor[ft]);
//        }
//    }
//
    double somatic_prior_tumor[DIGT_SGRID::PRESTRAND_SIZE];
    for(unsigned ft(0); ft<DIGT_SGRID::PRESTRAND_SIZE; ++ft)
        somatic_prior_tumor[ft] = 1.0/static_cast<double>(DIGT_SGRID::PRESTRAND_SIZE-1);

    for (unsigned ngt(0); ngt<DIGT_SIMPLE::SIZE; ++ngt)
    {
        // logP(Gn=ngt, Gt=tgt)
        double log_diploid_prior_prob = bare_lnprior_normal[ngt];  // logP(Gn)
        for (unsigned tgt(0); tgt<TWO_STATE_SOMATIC::SIZE; ++tgt) // 0: non-somatic, 1: somatic
        {
            double log_prior_prob = log_diploid_prior_prob + ((tgt == 0) ? lnmatch : lnmismatch);
            double max_log_sum = -INFINITY;
//            double log_sum[DDIGT_SGRID::SIZE];
            double log_sum[DDIGT_SGRID::PRESTRAND_SIZE];

            for (unsigned ft(0); ft<DIGT_SGRID::PRESTRAND_SIZE; ++ft)
            {
                for (unsigned fn(0); fn<DIGT_SGRID::PRESTRAND_SIZE; ++fn)
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
                           lprob_f_given_g = -INFINITY;
                       }
                    }
                    else    // somatic
                    {
                       if(fn == ft)
                       {
                           lprob_f_given_g = -INFINITY;
                       }
                       else
                       {
                           if(!is_normal_contaminated || ngt != DIGT_SIMPLE::REF)
                           {
                               lprob_f_given_g = log_error_mod + ((fn == ngt) ? ln_cse_rate : ln_se_rate+log_error_mod);
                           }
                           else
                           {
                               // fn should be smaller than ft
                               if(get_fraction_from_index_snv(fn) >= 0.15*get_fraction_from_index_snv(ft))
                                   lprob_f_given_g = -INFINITY;
                               else
                               {
                                   lprob_f_given_g = std::log(somatic_prior_normal[fn]) + std::log(somatic_prior_tumor[ft]);
                               }
                           }
                       }
                    }

                    const unsigned dgt(DDIGT_SGRID::get_state(fn, ft));

                    double sum = lprob_f_given_g + normal_lhood[fn] + tumor_lhood[ft];
                    log_sum[dgt] = sum;

//                    printf("%d\t%d\t%lf\n", fn, ft, log_prior_prob + lprob_f_given_g);
                    if(sum > max_log_sum) max_log_sum = sum;
                }
            }

            // Currently, the single-strand noise states are ignored.
            // Now consider the single-strand noise states. note that these states
            // are unique in that we don't look for mixtures of somatic
            // variation with these noise states, b/c single-strand
            // observations can almost exclusively be ruled out as noise:
            //
//            for (unsigned f(DIGT_SGRID::PRESTRAND_SIZE); f<DIGT_SGRID::SIZE; ++f)
//            {
//                const unsigned dgt(DDIGT_SGRID::get_state(f,f));
////                if(tgt == DIGT_SIMPLE::HET)
////                {
////                    double lprob_f_given_g = bare_lnprior_normal[ngt] + lnmatch + ln_se_rate+log_error_mod;
////                    double sum = normal_lhood[f]+tumor_lhood[f]+lprob_f_given_g;
////                    log_sum[dgt] = sum;
////                    if(sum > max_log_sum) max_log_sum = sum;
////                }
////                else
////                {
////                    log_sum[dgt] = -INFINITY;
////                }
//                log_sum[dgt] = -INFINITY;
//            }

            // calculate log(exp(log_sum[0])+exp(log_sum[1])+...)
            double sum = 0.0;
//            for (int i(0); i<DDIGT_SGRID::SIZE; ++i)
            for (int i(0); i<DDIGT_SGRID::PRESTRAND_SIZE; ++i)
            {
                sum += std::exp(log_sum[i] - max_log_sum);
            }

            log_post_prob[ngt][tgt] = log_prior_prob + max_log_sum + std::log(sum);

            if(log_post_prob[ngt][tgt] > max_log_prob)
            {
                max_log_prob = log_post_prob[ngt][tgt];
            }
        }
    }

    // Calculate posterior probabilities
    double sum_prob = 0.0;
    for (unsigned ngt(0); ngt<DIGT_SIMPLE::SIZE; ++ngt)
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

    double post_prob[DIGT_SIMPLE::SIZE][TWO_STATE_SOMATIC::SIZE];
    for (unsigned ngt(0); ngt<DIGT_SIMPLE::SIZE; ++ngt)
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
//            printf("%d\t%d\t%lf\n", ngt, tgt, post_prob[ngt][tgt]);
        }

        double err_som_and_ngt = 1.0 - som_prob_given_ngt;
        if (err_som_and_ngt < min_not_somfrom_sum)
        {
            min_not_somfrom_sum=err_som_and_ngt;
            rs.snv_from_ntype_qphred=error_prob_to_qphred(err_som_and_ngt);
            rs.ntype=ngt;
        }
    }
    rs.snv_qphred=error_prob_to_qphred(nonsom_prob);

    if ((! (is_forced_output || isComputeNonSomatic)) && (0==rs.snv_qphred)) return;

    // TODO: Calculate QSS_NT assuming polymorphic prior

    //
    // add new somatic gVCF value -- note this is an expanded definition of 'non-somatic' beyond just f_N == f_T
    // process regular tumor/normal lhood, but:
    // (1) use uniform probability for {somatic,non-somatic} states
    // (2) simply computation to remove strand-specific logic
    // (3) ignore normal genotype
    rs.nonsomatic_qphred=rs.snv_qphred; // TODO: not implemented yet

    static const bool is_compute_sb(true);
    if (is_compute_sb)
    {
        // get ratio of strand bias vs. non-strand-bias version of max_gt, if max_gt does not correspond to a het state, then
        // set sb to 0
        const blt_float_t symm_lhood(*std::max_element(tumor_lhood+N_BASE, tumor_lhood+DIGT_SGRID::PRESTRAND_SIZE));
        const blt_float_t strand_lhood(*std::max_element(tumor_lhood+DIGT_SGRID::PRESTRAND_SIZE, tumor_lhood+DIGT_SGRID::SIZE));
        rs.strandBias = std::max(0.f,(strand_lhood - symm_lhood));
    }
}

void
somatic_snv_caller_strand_grid::
position_somatic_snv_call(
    const extended_pos_info& normal_epi,
    const extended_pos_info& tumor_epi,
    const extended_pos_info* normal_epi_t2_ptr,
    const extended_pos_info* tumor_epi_t2_ptr,
    const bool isComputeNonSomatic,
    somatic_snv_genotype_grid& sgt) const
{
    {
        const snp_pos_info& normal_pi(normal_epi.pi);
        const snp_pos_info& tumor_pi(tumor_epi.pi);

        if (normal_pi.get_ref_base()=='N')
        {
            sgt.is_forced_output=false;
            return;
        }
        sgt.ref_gt=base_to_id(normal_pi.get_ref_base());

        // check that a non-reference call meeting quality criteria even
        // exists:
        if (! (sgt.is_forced_output || isComputeNonSomatic))
        {
            if (is_spi_allref(normal_pi,sgt.ref_gt) && is_spi_allref(tumor_pi,sgt.ref_gt)) return;
        }
    }

    // strawman model treats normal and tumor as independent, so
    // calculate separate lhoods:
    blt_float_t normal_lhood[DIGT_SGRID::SIZE];
    blt_float_t tumor_lhood[DIGT_SGRID::SIZE];

    const bool is_tier2(NULL != normal_epi_t2_ptr);

    static const unsigned n_tier(2);
    result_set tier_rs[n_tier];
    for (unsigned i(0); i<n_tier; ++i)
    {
        const bool is_include_tier2(i==1);
        if (is_include_tier2)
        {
            if (! is_tier2) continue;
            if (tier_rs[0].snv_qphred==0)
            {
                tier_rs[1] = tier_rs[0];
                continue;
            }
        }

        // get likelihood of each genotype (REF, HOM, HET)
        const extended_pos_info& nepi(is_include_tier2 ? *normal_epi_t2_ptr : normal_epi );
        const extended_pos_info& tepi(is_include_tier2 ? *tumor_epi_t2_ptr : tumor_epi );

        // determine alt base (experimental)
        // later, this must be included in pi
        int alt_count[4] = {};
        for (const base_call& tbc : tepi.pi.calls)
        {
            const uint8_t obs_id(tbc.base_id);
            if(obs_id == sgt.ref_gt || obs_id >= 4) continue;
            ++alt_count[obs_id];
        }
        unsigned alt_id = sgt.ref_gt;
        int max_count = 0;
        for (unsigned idx(0); idx<4; ++idx)
        {
            if (alt_count[idx] > max_count)
            {
                max_count = alt_count[idx];
                alt_id = idx;
            }
        }

        get_diploid_gt_lhood_cached_simple(nepi.pi, sgt.ref_gt, normal_lhood);
        get_diploid_gt_lhood_cached_simple(tepi.pi, sgt.ref_gt, tumor_lhood);

        // get likelihood of non-canonical frequencies (0.05, 0.1, ..., 0.45, 0.55, ..., 0.95)
        get_diploid_het_grid_lhood_cached(nepi.pi, sgt.ref_gt, DIGT_SGRID::HET_RES, normal_lhood+DIGT_SIMPLE::SIZE);
        get_diploid_het_grid_lhood_cached(tepi.pi, sgt.ref_gt, DIGT_SGRID::HET_RES, tumor_lhood+DIGT_SIMPLE::SIZE);

        // get likelihood of strand states (0.05, ..., 0.45)
        get_diploid_strand_grid_lhood_spi(nepi.pi,sgt.ref_gt,normal_lhood+DIGT_SGRID::PRESTRAND_SIZE);
        get_diploid_strand_grid_lhood_spi(tepi.pi,sgt.ref_gt,tumor_lhood+DIGT_SGRID::PRESTRAND_SIZE);

        // genomic site results:
        calculate_result_set_grid(isComputeNonSomatic,
                                  normal_lhood,
                                  tumor_lhood,
                                  _bare_lnprior,
                                  _ln_nostrand_sse_rate,
                                  _ln_csse_rate,
                                  _ln_som_match,_ln_som_mismatch,
                                  sgt.is_forced_output,
                                  tier_rs[i]);
        tier_rs[i].max_gt = alt_id;


#if 0
#ifdef ENABLE_POLY
        // polymorphic site results:
        assert(0); // still needs to be adapted for 2-tier system:
        calculate_result_set(normal_lhood,tumor_lhood,
                             lnprior_polymorphic(sgt.ref_gt),sgt.ref_gt,sgt.poly);
#else
        sgt.poly.snv_qphred = 0;
#endif
#endif

#ifdef SOMATIC_DEBUG
        if ((i==0) && ((tier_rs[i].snv_qphred > 0) || isComputeNonSomatic))
        {
            const somatic_snv_caller_strand_grid::prior_set& pset(get_prior_set(sgt.ref_gt));
#ifdef STANDARD
            const blt_float_t lnmatch(_ln_som_match);
            const blt_float_t lnmismatch(_ln_som_mismatch);

            const unsigned state_size(DDIGT_SGRID::SIZE);
#else
            const unsigned state_size(DDIGT_SGRID::PRESTRAND_SIZE);
#endif

            log_os << "DUMP ON\n";
            log_os << "tier1_qphred_snv: " << tier_rs[0].snv_qphred << "\n";
            log_os << "tier1_qphred_nonsomatic: " << tier_rs[0].nonsomatic_qphred << "\n";

            // instead of dumping the entire distribution, we sort the lhood,prior,and prob to print out the N top values of each:
            std::vector<double> lhood(state_size,0);
            std::vector<double> prior(state_size,0);
            std::vector<double> post(state_size,0);

            // first get raw lhood:
            //
            for (unsigned ngt(0); ngt<DIGT_SGRID::PRESTRAND_SIZE; ++ngt)
            {
                for (unsigned tgt(0); tgt<DIGT_SGRID::PRESTRAND_SIZE; ++tgt)
                {
                    const unsigned dgt(DDIGT_SGRID::get_state(ngt,tgt));
                    // unorm takes the role of the normal prior for the somatic case:
                    //            static const blt_float_t unorm(std::log(static_cast<blt_float_t>(DIGT_SGRID::PRESTRAND_SIZE)));

                    // switch between standard and gvcf info:
#ifdef STANDARD
                    //blt_float_t prior;
                    //if(tgt==ngt) { prior=pset.normal[ngt]+lnmatch; }
                    //else         { prior=pset.somatic_marginal[ngt]+lnmismatch; }
                    blt_float_t pr;
                    if (tgt==ngt)
                    {
                        pr=pset.normal[ngt]+lnmatch;
                    }
                    else
                    {
                        pr=pset.somatic_marginal[ngt]+lnmismatch;
                    }
                    prior[dgt] = pr;
#else
                    prior[dgt] = nonsomatic_gvcf_prior(pset,ngt,tgt);
#endif
                    lhood[dgt] = normal_lhood[ngt]+tumor_lhood[tgt];
                    post[dgt] = lhood[dgt] + prior[dgt];
                }
            }

#ifdef STANDARD
            for (unsigned gt(DIGT_SGRID::PRESTRAND_SIZE); gt<DIGT_SGRID::SIZE; ++gt)
            {
                const unsigned dgt(DDIGT_SGRID::get_state(gt,gt));
                lhood[dgt] = normal_lhood[gt]+tumor_lhood[gt];
                prior[dgt] = pset.normal[gt]+lnmatch;
                post[dgt] = lhood[dgt] + prior[dgt];
            }
#endif

            std::vector<double> lhood2(lhood);
            sort_n_dump("lhood_prior",lhood,prior,sgt.ref_gt);
            sort_n_dump("post_lhood",post,lhood2,sgt.ref_gt);

            log_os << "DUMP OFF\n";
        }
#endif

    }

    if (! (sgt.is_forced_output || isComputeNonSomatic))
    {
        if ((tier_rs[0].snv_qphred==0) ||
            (is_tier2 && (tier_rs[1].snv_qphred==0))) return;
    }

    sgt.snv_tier=0;
    sgt.snv_from_ntype_tier=0;
    if (is_tier2)
    {
        if (tier_rs[0].snv_qphred > tier_rs[1].snv_qphred)
        {
            sgt.snv_tier=1;
        }

        if (tier_rs[0].snv_from_ntype_qphred > tier_rs[1].snv_from_ntype_qphred)
        {
            sgt.snv_from_ntype_tier=1;
        }
    }

    sgt.rs=tier_rs[sgt.snv_from_ntype_tier];

    if (is_tier2 && (tier_rs[0].ntype != tier_rs[1].ntype))
    {
        // catch NTYPE conflict states:
        sgt.rs.ntype = NTYPE::CONFLICT;
        sgt.rs.snv_from_ntype_qphred = 0;
    }
    else
    {
        // classify NTYPE:
        //

        // convert diploid genotype into more limited ntype set:
        //
        if       (sgt.rs.ntype==DIGT_SIMPLE::REF)
        {
            sgt.rs.ntype=NTYPE::REF;
        }
        else if (sgt.rs.ntype==DIGT_SIMPLE::HOM)
        {
            sgt.rs.ntype=NTYPE::HOM;
        }
        else
        {
            sgt.rs.ntype=NTYPE::HET;
        }
    }

    sgt.rs.snv_qphred = tier_rs[sgt.snv_tier].snv_qphred;

    /// somatic gVCF, always use tier1 to keep things simple:
    sgt.rs.nonsomatic_qphred = tier_rs[0].nonsomatic_qphred;
}
