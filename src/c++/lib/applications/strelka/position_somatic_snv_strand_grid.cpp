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

inline int get_prior_index(unsigned ngt, unsigned tgt, unsigned ft, unsigned fn)
{
    return ngt*TWO_STATE_SOMATIC::SIZE*DIGT_SGRID::PRESTRAND_SIZE*DIGT_SGRID::PRESTRAND_SIZE
            + tgt*DIGT_SGRID::PRESTRAND_SIZE*DIGT_SGRID::PRESTRAND_SIZE
            + ft*DIGT_SGRID::PRESTRAND_SIZE
            + fn;
}

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
get_prior(
        const blt_float_t ln_se_rate,   // ln (somatic error rate)
        const blt_float_t ln_cse_rate,  // ln (1 - somatic_error_rate)
        std::vector<blt_float_t>& ln_freq_given_somatic
        )
{
    bool is_normal_contaminated = true;
    // Calculate posterior probabilities

    blt_float_t log_error_mod = -std::log(static_cast<double>(DIGT_SGRID::PRESTRAND_SIZE-1));

    double somatic_prior_normal[DIGT_SGRID::PRESTRAND_SIZE] = {};
    somatic_prior_normal[DIGT_SIMPLE::REF] = 0.5;
    somatic_prior_normal[DIGT_SGRID::PRESTRAND_SIZE - 1] = 0.5;

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
        for (unsigned tgt(0); tgt<TWO_STATE_SOMATIC::SIZE; ++tgt) // 0: non-somatic, 1: somatic
        {
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

                    int index = get_prior_index(ngt, tgt, ft, fn);
                    ln_freq_given_somatic[index] = lprob_f_given_g;
                }
            }
        }
    }
}

somatic_snv_caller_strand_grid::
somatic_snv_caller_strand_grid(const strelka_options& opt,
                               const pprob_digt_caller& pd_caller)
: _ln_freq_given_somatic(DIGT_SIMPLE::SIZE*TWO_STATE_SOMATIC::SIZE*DIGT_SGRID::PRESTRAND_SIZE*DIGT_SGRID::PRESTRAND_SIZE),
  _ln_som_match(log1p_switch(-opt.somatic_snv_rate)),
  _ln_som_mismatch(std::log(opt.somatic_snv_rate))
{
    _bare_lnprior = pd_caller.lnprior_genomic();

    const blt_float_t strand_sse_rate(opt.shared_site_error_rate*opt.shared_site_error_strand_bias_fraction);
    const blt_float_t nostrand_sse_rate(opt.shared_site_error_rate-strand_sse_rate);

    blt_float_t ln_csse_rate = log1p_switch(-opt.shared_site_error_rate);
    blt_float_t ln_nostrand_sse_rate = std::log(nostrand_sse_rate);

    std::fill(_ln_freq_given_somatic.begin(),_ln_freq_given_somatic.end(),0);

    get_prior(
            ln_nostrand_sse_rate,   // ln (somatic error rate)
            ln_csse_rate,  // ln (1 - somatic_error_rate)
            _ln_freq_given_somatic);
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

static
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
        )
{
    double log_post_prob[DIGT_SIMPLE::SIZE][TWO_STATE_SOMATIC::SIZE];
    double max_log_prob = -INFINITY;


    for (unsigned ngt(0); ngt<DIGT_SIMPLE::SIZE; ++ngt)
    {
        // logP(Gn=ngt, Gt=tgt)
        double log_diploid_prior_prob = bare_lnprior_normal[ngt];  // logP(Gn)
        for (unsigned tgt(0); tgt<TWO_STATE_SOMATIC::SIZE; ++tgt) // 0: non-somatic, 1: somatic
        {
            double log_prior_prob = log_diploid_prior_prob + ((tgt == 0) ? lnmatch : lnmismatch);
            double max_log_sum = -INFINITY;
            double log_sum[DDIGT_SGRID::PRESTRAND_SIZE];

            for (unsigned ft(0); ft<DIGT_SGRID::PRESTRAND_SIZE; ++ft)
            {
                for (unsigned fn(0); fn<DIGT_SGRID::PRESTRAND_SIZE; ++fn)
                {
                    const int prior_index = get_prior_index(ngt, tgt, ft, fn);

                    const unsigned dgt(DDIGT_SGRID::get_state(fn, ft));

                    double sum = ln_freq_given_somatic[prior_index] + normal_lhood[fn] + tumor_lhood[ft];
                    log_sum[dgt] = sum;

//                    printf("%d\t%d\t%lf\n", fn, ft, log_prior_prob + lprob_f_given_g);
                    if(sum > max_log_sum) max_log_sum = sum;
                }
            }

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
        const blt_float_t symm_lhood(*std::max_element(tumor_lhood+DIGT_SIMPLE::SIZE, tumor_lhood+DIGT_SGRID::PRESTRAND_SIZE));
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

        // determine alt id
        // TODO: consider to include it in pi
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
//        get_diploid_strand_grid_lhood_spi(nepi.pi,sgt.ref_gt,normal_lhood+DIGT_SGRID::PRESTRAND_SIZE);
        get_diploid_strand_grid_lhood_spi(tepi.pi,sgt.ref_gt,tumor_lhood+DIGT_SGRID::PRESTRAND_SIZE);

        // genomic site results:
        calculate_result_set_grid(isComputeNonSomatic,
                                  normal_lhood,
                                  tumor_lhood,
                                  _ln_freq_given_somatic,
                                  _bare_lnprior,
                                  _ln_som_match,_ln_som_mismatch,
                                  sgt.is_forced_output,
                                  tier_rs[i]);
        tier_rs[i].max_gt = alt_id;
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
