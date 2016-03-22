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
#include "position_somatic_snv_strand_grid_lhood_cached.hh"
#include "qscore_calculator.hh"
#include "somatic_result_set.hh"
#include "somatic_call_shared.hh"
#include "blt_common/snp_util.hh"
#include "blt_util/log.hh"
#include "blt_util/math_util.hh"
#include "blt_util/prob_util.hh"
#include "blt_util/seq_util.hh"

#include "strelka_common/het_ratio_cache.hh"

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

somatic_snv_caller_strand_grid::
somatic_snv_caller_strand_grid(const strelka_options& opt)
: _ln_somatic_prior(SOMATIC_DIGT::SIZE*TWO_STATE_SOMATIC::SIZE*DIGT_GRID::PRESTRAND_SIZE*DIGT_GRID::PRESTRAND_SIZE),
  _ln_som_match(log1p_switch(-opt.somatic_snv_rate)),
  _ln_som_mismatch(std::log(opt.somatic_snv_rate))
{
    calculate_bare_lnprior(opt.bsnp_diploid_theta, _bare_lnprior);

    const blt_float_t strand_sse_rate(opt.shared_site_error_rate*opt.shared_site_error_strand_bias_fraction);
    const blt_float_t nostrand_sse_rate(opt.shared_site_error_rate-strand_sse_rate);

    blt_float_t ln_csse_rate = log1p_switch(-opt.shared_site_error_rate);
    blt_float_t ln_sse_rate = std::log(nostrand_sse_rate);

    std::fill(_ln_somatic_prior.begin(),_ln_somatic_prior.end(),0);

    set_prior(
            opt.ssnv_freq_ratio,
            ln_sse_rate,   // ln (shared_error_rate)
            ln_csse_rate,  // ln (1 - shared_error_rate)
            _ln_somatic_prior);
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
    blt_float_t lhood_fwd = 0; // "on-strand" is fwd
    blt_float_t lhood_rev = 0; // "on-strand" is rev

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
            lhood_fwd += val_fwd;
            lhood_rev += val_rev;
        }
        else
        {
            const blt_float_t val_off_strand(bc.ln_error_prob()+ln_one_third);
            const blt_float_t val_fwd(bc.is_fwd_strand ? cv.val[1] : val_off_strand);
            const blt_float_t val_rev(bc.is_fwd_strand ? val_off_strand : cv.val[1]);

            lhood_fwd += val_fwd;
            lhood_rev += val_rev;
        }
    }

    *lhood = log_sum(lhood_fwd,lhood_rev)+ln_one_half;
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
    for (unsigned gt(0); gt<(DIGT_GRID::STRAND_STATE_SIZE); ++gt) lhood[gt] = 0.;

    static const blt_float_t ratio_increment(0.5/static_cast<blt_float_t>(DIGT_GRID::HET_RES+1));
    for (unsigned i(0); i<DIGT_GRID::HET_RES; ++i)
    {
        const blt_float_t het_ratio((i+1)*ratio_increment);
        get_strand_ratio_lhood_spi(pi,ref_gt,het_ratio,i,hrcache,
                                   lhood+i);
    }
}

/// expanded definition of 'nonsomatic' for the purpose of providing somatic gVCF output:
static
float
gvcf_nonsomatic_gvcf_prior(
    const unsigned fn,
    const unsigned ft)
{
    if (fn == ft)
    {
        static const float lone(std::log(1.f));
        return lone;
    }
    else if (fn == SOMATIC_DIGT::REF || fn == SOMATIC_DIGT::HET)
    {
        static const float lhalf(std::log(0.5f));
        return lhalf;
    }
    else
    {
        static const float lzero(-std::numeric_limits<float>::infinity());
        return lzero;
    }
}

static
void
calculate_result_set_grid(
        const bool isComputeNonSomatic,
        const std::vector<blt_float_t>& ln_somatic_prior,
        const blt_float_t* normal_lhood,
        const blt_float_t* tumor_lhood,
        const blt_float_t* bare_lnprior_normal,
        const blt_float_t lnmatch,
        const blt_float_t lnmismatch,
        const bool is_forced_output,
        snv_result_set& rs
        )
{

    calculate_result_set_grid(
            normal_lhood,
            tumor_lhood,
            ln_somatic_prior,
            bare_lnprior_normal,
            lnmatch,
            lnmismatch,
            rs);

    if ((! (is_forced_output || isComputeNonSomatic)) && (0==rs.qphred)) return;

    // add new somatic gVCF value -- note this is an expanded definition of 'non-somatic' beyond just f_N == f_T
    if (isComputeNonSomatic)
    {
        // process regular tumor/normal lhood, but:
        // (1) use uniform probability for {somatic,non-somatic} states
        // (2) simply computation to remove strand-specific logic
        // (3) ignore normal genotype
        //
        std::vector<double> pprob(DDIGT_GRID::SIZE);
        double sum_prob = 0.0;
        for (unsigned fn(0); fn<DIGT_GRID::PRESTRAND_SIZE; ++fn)
        {
            for (unsigned ft(0); ft<DIGT_GRID::PRESTRAND_SIZE; ++ft)
            {
                const unsigned dgt(DDIGT_GRID::get_state(fn,ft));
                double prob = normal_lhood[fn]+tumor_lhood[ft]+gvcf_nonsomatic_gvcf_prior(fn,ft);
                pprob[dgt] = prob;
                sum_prob += prob;
            }
        }

        double sgvcf_nonsomatic_sum(0);
        for (unsigned f(0); f<DIGT_GRID::PRESTRAND_SIZE; ++f)
        {
            const unsigned dgt(DDIGT_GRID::get_state(f,f));
            sgvcf_nonsomatic_sum += pprob[dgt] / sum_prob;
        }

        rs.nonsomatic_qphred=error_prob_to_qphred(1.-sgvcf_nonsomatic_sum);
    }

    static const bool is_compute_sb(true);
    if (is_compute_sb)
    {
        // get ratio of strand bias vs. non-strand-bias version of max_gt, if max_gt does not correspond to a het state, then
        // set sb to 0
        const blt_float_t symm_lhood(*std::max_element(tumor_lhood+SOMATIC_DIGT::SIZE, tumor_lhood+DIGT_GRID::PRESTRAND_SIZE));
        const blt_float_t strand_lhood(*std::max_element(tumor_lhood+DIGT_GRID::PRESTRAND_SIZE, tumor_lhood+DIGT_GRID::SIZE));
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
    blt_float_t normal_lhood[DIGT_GRID::SIZE];
    blt_float_t tumor_lhood[DIGT_GRID::SIZE];

    const bool is_tier2(NULL != normal_epi_t2_ptr);

    static const unsigned n_tier(2);
    snv_result_set tier_rs[n_tier];
    for (unsigned i(0); i<n_tier; ++i)
    {
        const bool is_include_tier2(i==1);
        if (is_include_tier2)
        {
            if (! is_tier2) continue;
            if (tier_rs[0].qphred==0)
            {
                tier_rs[1] = tier_rs[0];
                continue;
            }
        }

        // get likelihood of each genotype (REF, HOM, HET)
        const extended_pos_info& nepi(is_include_tier2 ? *normal_epi_t2_ptr : normal_epi );
        const extended_pos_info& tepi(is_include_tier2 ? *tumor_epi_t2_ptr : tumor_epi );

        get_diploid_gt_lhood_cached_simple(nepi.pi, sgt.ref_gt, normal_lhood);
        get_diploid_gt_lhood_cached_simple(tepi.pi, sgt.ref_gt, tumor_lhood);

        // get likelihood of non-canonical frequencies (0.05, 0.1, ..., 0.45, 0.55, ..., 0.95)
        get_diploid_het_grid_lhood_cached(nepi.pi, sgt.ref_gt, DIGT_GRID::HET_RES, normal_lhood+SOMATIC_DIGT::SIZE);
        get_diploid_het_grid_lhood_cached(tepi.pi, sgt.ref_gt, DIGT_GRID::HET_RES, tumor_lhood+SOMATIC_DIGT::SIZE);

        // get likelihood of strand states (0.05, ..., 0.45)
//        get_diploid_strand_grid_lhood_spi(nepi.pi,sgt.ref_gt,normal_lhood+DIGT_GRID::PRESTRAND_SIZE);
        get_diploid_strand_grid_lhood_spi(tepi.pi,sgt.ref_gt,tumor_lhood+DIGT_GRID::PRESTRAND_SIZE);

        // genomic site results:
        calculate_result_set_grid(isComputeNonSomatic,
                                  _ln_somatic_prior,
                                  normal_lhood,
                                  tumor_lhood,
                                  _bare_lnprior,
                                  _ln_som_match,_ln_som_mismatch,
                                  sgt.is_forced_output,
                                  tier_rs[i]);
    }

    if (! (sgt.is_forced_output || isComputeNonSomatic))
    {
        if ((tier_rs[0].qphred==0) ||
            (is_tier2 && (tier_rs[1].qphred==0))) return;
    }

    sgt.snv_tier=0;
    sgt.snv_from_ntype_tier=0;
    if (is_tier2)
    {
        if (tier_rs[0].qphred > tier_rs[1].qphred)
        {
            sgt.snv_tier=1;
        }

        if (tier_rs[0].from_ntype_qphred > tier_rs[1].from_ntype_qphred)
        {
            sgt.snv_from_ntype_tier=1;
        }
    }

    sgt.rs=tier_rs[sgt.snv_from_ntype_tier];

    if (is_tier2 && (tier_rs[0].ntype != tier_rs[1].ntype))
    {
        // catch NTYPE conflict states:
        sgt.rs.ntype = NTYPE::CONFLICT;
        sgt.rs.from_ntype_qphred = 0;
    }
    else
    {
        // classify NTYPE:
        //

        // convert diploid genotype into more limited ntype set:
        //
        if       (sgt.rs.ntype==SOMATIC_DIGT::REF)
        {
            sgt.rs.ntype=NTYPE::REF;
        }
        else if (sgt.rs.ntype==SOMATIC_DIGT::HOM)
        {
            sgt.rs.ntype=NTYPE::HOM;
        }
        else
        {
            sgt.rs.ntype=NTYPE::HET;
        }
    }

    sgt.rs.qphred = tier_rs[sgt.snv_tier].qphred;

    /// somatic gVCF, always use tier1 to keep things simple:
    sgt.rs.nonsomatic_qphred = tier_rs[0].nonsomatic_qphred;

    sgt.rs.normal_alt_id = normal_epi.pi.get_most_frequent_alt_id(sgt.ref_gt);
    sgt.rs.tumor_alt_id = tumor_epi.pi.get_most_frequent_alt_id(sgt.ref_gt);
}
