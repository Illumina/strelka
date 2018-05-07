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
/// \author Sangtae Kim
///

#include "position_somatic_snv_strand_grid_lhood_cached.hh"
#include "strelka_digt_states.hh"
#include "strelka_common/het_ratio_cache.hh"

#include "blt_util/digt.hh"
#include "blt_util/logSumUtil.hh"

#include <cmath>


static const blt_float_t one_third(1./3.);
static const blt_float_t ln_one_third(std::log(one_third));
static const blt_float_t one_half(1./2.);
static const blt_float_t ln_one_half(std::log(one_half));


void
get_diploid_gt_lhood_cached_simple(
    const snp_pos_info& pi,
    const unsigned ref_gt,
    blt_float_t* const lhood)
{
    // ! not thread-safe !
    static het_ratio_cache<3> hrcache;

    // get likelihood of each genotype
    for (unsigned gt(0); gt<SOMATIC_DIGT::SIZE; ++gt) lhood[gt] = 0.;

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

        if (obs_id == ref_gt)
        {
            lhood[SOMATIC_DIGT::REF] += cv.val[2];
            lhood[SOMATIC_DIGT::HET] += cv.val[1];
            lhood[SOMATIC_DIGT::HOM] += cv.val[0];
        }
        else
        {
            lhood[SOMATIC_DIGT::REF] += cv.val[0];
            lhood[SOMATIC_DIGT::HET] += cv.val[1];
            lhood[SOMATIC_DIGT::HOM] += cv.val[2];
        }
    }
}

static
void
get_high_low_het_ratio_lhood_cached(
    const snp_pos_info& pi,
    const unsigned ref_gt,
    const blt_float_t het_ratio,
    const unsigned het_ratio_index,
    het_ratio_cache<2>& hrcache,
    blt_float_t* lhood_high,
    blt_float_t* lhood_low)
{
    const blt_float_t chet_ratio(1.-het_ratio);

    const unsigned n_calls(pi.calls.size());

    for (unsigned i(0); i<n_calls; ++i)
    {
        const base_call& bc(pi.calls[i]);

        std::pair<bool,cache_val<2>*> ret(hrcache.get_val(bc.get_qscore(),het_ratio_index));
        cache_val<2>& cv(*ret.second);
        if (! ret.first)
        {
            const blt_float_t eprob(bc.error_prob());
            const blt_float_t ceprob(1-eprob);

            // precalculate the result for expect values of het_ratio and chet_ratio
            cv.val[0] = std::log((ceprob)*het_ratio+((eprob)*one_third)*chet_ratio);    // mismatch for lhood_low, match for lhood_high
            cv.val[1] = std::log((ceprob)*chet_ratio+((eprob)*one_third)*het_ratio);    // match for lhood_low, mismatch for lhood_high
        }

        const uint8_t obs_id(bc.base_id);

        if (obs_id == ref_gt)   // match
        {
            *lhood_high += cv.val[0];
            *lhood_low += cv.val[1];
        }
        else                    // mismatch
        {
            *lhood_high += cv.val[1];
            *lhood_low += cv.val[0];
        }
    }
}

void
get_diploid_het_grid_lhood_cached(
    const snp_pos_info& pi,
    const unsigned ref_gt,
    const unsigned hetResolution,
    blt_float_t* const lhood)
{
    // ! not thread-safe !
    static het_ratio_cache<2> hrcache;

    // get likelihood of each genotype
    const unsigned totalHetRatios(hetResolution*2);
    for (unsigned gt(0); gt<totalHetRatios; ++gt) lhood[gt] = 0.;

//    blt_float_t* lhood_off=lhood-N_BASE;

    for (unsigned hetIndex(0); hetIndex<hetResolution; ++hetIndex)
    {
        const blt_float_t het_ratio((hetIndex+1)*DIGT_GRID::RATIO_INCREMENT);
        get_high_low_het_ratio_lhood_cached(pi,ref_gt, het_ratio,hetIndex,hrcache,
                                            lhood+(totalHetRatios-(hetIndex+1)),
                                            lhood+hetIndex);
    }
}

// calculate probability of strand-specific noise
//
// accelerated version with no hyrax q-val mods:
//
// the ratio key can be used as a proxy for the het ratio to look up cached results:
//
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

    *lhood = getLogSum(lhood_fwd,lhood_rev)+ln_one_half;
}

