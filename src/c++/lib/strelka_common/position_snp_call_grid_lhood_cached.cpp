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

/// \author Chris Saunders
///

#include "position_snp_call_grid_lhood_cached.hh"

#include "blt_util/digt.hh"
#include "blt_util/logSumUtil.hh"

#include <cmath>


static const blt_float_t one_third(1./3.);
static const blt_float_t ln_one_third(std::log(one_third));
static const blt_float_t one_half(1./2.);
static const blt_float_t ln_one_half(std::log(one_half));



// accelerated version with no hyrax q-val mods:
//
// the ratio key can be used as a proxy for the het ratio to look up cached results:
//
static
void
get_high_low_het_ratio_lhood_cached(
    const snp_pos_info& pi,
    const blt_float_t het_ratio,
    const unsigned het_ratio_index,
    het_ratio_cache<3>& hrcache,
    blt_float_t* lhood_high,
    blt_float_t* lhood_low)
{
    const blt_float_t chet_ratio(1.-het_ratio);

    const unsigned n_calls(pi.calls.size());

    //    cache_val cv;
    static const uint8_t remap[3] = {0,2,1};

    for (unsigned i(0); i<n_calls; ++i)
    {
        const base_call& bc(pi.calls[i]);

        std::pair<bool,cache_val<3>*> ret(hrcache.get_val(bc.get_qscore(),het_ratio_index));
        cache_val<3>& cv(*ret.second);
        if (! ret.first)
        {
            const blt_float_t eprob(bc.error_prob());
            const blt_float_t ceprob(1-eprob);
            //const blt_float_t lne(bc.ln_error_prob());
            //const blt_float_t lnce(bc.ln_comp_error_prob());

            // precalculate the result for expect values of 0.0, het_ratio, chet_ratio, 1.0
            cv.val[0] = bc.ln_error_prob()+ln_one_third;
            cv.val[1] = std::log((ceprob)*het_ratio+((eprob)*one_third)*chet_ratio);
            cv.val[2] = std::log((ceprob)*chet_ratio+((eprob)*one_third)*het_ratio);
        }

        const uint8_t obs_id(bc.base_id);

        for (unsigned gt(N_BASE); gt<DIGT::SIZE; ++gt)
        {
            const unsigned key(DIGT::expect2_bias(obs_id,gt));
            lhood_high[gt] += cv.val[key];
            lhood_low[gt] += cv.val[remap[key]];
        }
    }
}



// accelerated version with no hyrax q-val mods:
//
static
void
increment_het_ratio_lhood_cached(
    const snp_pos_info& pi,
    const blt_float_t het_ratio,
    const unsigned het_ratio_index,
    het_ratio_cache<3>& hrcache,
    blt_float_t* all_het_lhood)
{
    // multiply probs of alternate ratios into local likelihoods, then
    // *add* them to the global tally (effectively this is the sum lhood of
    // many different heterozygous genotypes).
    //
    // in the gt_high genotype, the first allele (in lexicographical
    // order) is expected at het_ratio and the second allele is
    // expected at chet_ratio.  gt_low genotype is vice versa.
    //
    blt_float_t lhood_high[DIGT::SIZE];
    blt_float_t lhood_low[DIGT::SIZE];
    for (unsigned gt(0); gt<DIGT::SIZE; ++gt)
    {
        lhood_high[gt] = 0.;
        lhood_low[gt] = 0.;
    }
    get_high_low_het_ratio_lhood_cached(pi,het_ratio,het_ratio_index,hrcache,lhood_high,lhood_low);

    for (unsigned gt(0); gt<DIGT::SIZE; ++gt)
    {
        if (! DIGT::is_het(gt)) continue;
        all_het_lhood[gt] = getLogSum(all_het_lhood[gt], lhood_high[gt], lhood_low[gt]);
    }
}



void
get_diploid_gt_lhood_cached(
    const blt_options& opt,
    const snp_pos_info& pi,
    const bool useHetVariantFrequencyExtension,
    const blt_float_t hetVariantFrequencyExtension,
    blt_float_t* const lhood)
{
    // ! not thread-safe !
    static het_ratio_cache<3> hrcache;

    // get likelihood of each genotype
    for (unsigned gt(0); gt<DIGT::SIZE; ++gt) lhood[gt] = 0.;

    for (const base_call& bc : pi.calls)
    {
        std::pair<bool,cache_val<3>*> ret(hrcache.get_val(bc.get_qscore(),0));
        cache_val<3>& cv(*ret.second);
        if (! ret.first)
        {
            const blt_float_t eprob(bc.error_prob());
            const blt_float_t ceprob(1-eprob);
            const blt_float_t lne(bc.ln_error_prob());
            const blt_float_t lnce(bc.ln_comp_error_prob());

            // precalculate the result for expect values of 0.0, 0.5 & 1.0
            cv.val[0] = lne+ln_one_third;
            cv.val[1] = std::log((ceprob)+((eprob)*one_third))+ln_one_half;
            cv.val[2] = lnce;
        }

        const uint8_t obs_id(bc.base_id);
        for (unsigned gt(0); gt<DIGT::SIZE; ++gt)
        {
            lhood[gt] += cv.val[DIGT::expect2(obs_id,gt)];
        }
    }

    // het bias here refers to an expanded frequency range for the
    // heterozygous state, referred to as the myrax snp calling with in
    // single-sample analysis and not currently used for somatic
    // calls (as of strelka proto3/4)
    //
    if (useHetVariantFrequencyExtension)
    {
        // ! not thread-safe !
        static het_ratio_cache<3> hrcache_bias;

        // loop is currently setup to assume a uniform het ratio subgenotype prior
        const unsigned n_bias_steps(1+static_cast<unsigned>(hetVariantFrequencyExtension/opt.maxHetVariantFrequencyIncrement));
        const blt_float_t ratio_increment(hetVariantFrequencyExtension/static_cast<blt_float_t>(n_bias_steps));
        for (unsigned i(0); i<n_bias_steps; ++i)
        {
            const blt_float_t het_ratio(0.5+(i+1)*ratio_increment);
            increment_het_ratio_lhood_cached(pi,het_ratio,i,hrcache_bias,lhood);
        }

        const unsigned n_het_subgt(1+2*n_bias_steps);
        const blt_float_t subgt_log_prior(std::log(1./static_cast<blt_float_t>(n_het_subgt)));

        for (unsigned gt(0); gt<DIGT::SIZE; ++gt)
        {
            if (! DIGT::is_het(gt)) continue;
            lhood[gt] += subgt_log_prior;
        }
    }
}



// fill in noise portions of the likelihood distro for non-strand
// noise
//
void
get_diploid_het_grid_lhood_cached(
    const snp_pos_info& pi,
    const unsigned hetResolution,
    blt_float_t* const lhood)
{
    // ! not thread-safe !
    static het_ratio_cache<3> hrcache;

    // get likelihood of each genotype
    const unsigned totalHetRatios(hetResolution*2);
    for (unsigned gt(0); gt<(totalHetRatios*DIGT::HET_SIZE); ++gt) lhood[gt] = 0.;

    blt_float_t* lhood_off=lhood-N_BASE;

    const blt_float_t ratio_increment(0.5/static_cast<blt_float_t>(hetResolution+1));
    for (unsigned hetIndex(0); hetIndex<hetResolution; ++hetIndex)
    {
        const blt_float_t het_ratio((hetIndex+1)*ratio_increment);
        get_high_low_het_ratio_lhood_cached(pi,het_ratio,hetIndex,hrcache,
                                            lhood_off+(hetIndex*DIGT::HET_SIZE),
                                            lhood_off+((totalHetRatios-(hetIndex+1))*DIGT::HET_SIZE));
    }
}
