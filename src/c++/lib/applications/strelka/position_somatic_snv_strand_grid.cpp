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

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <map>

//#define SOMATIC_DEBUG

#ifdef SOMATIC_DEBUG
#include <iostream>
#include <iomanip>
#endif


constexpr blt_float_t one_third(1./3.);
static const blt_float_t ln_one_third(std::log(one_third));
constexpr blt_float_t one_half(1./2.);
static const blt_float_t ln_one_half(std::log(one_half));



// compute just the non-strand-bias portion of the normal marginal
// prior given p(signal), p(no-strand noise), p(strand-bias noise)
//
static
void
get_nostrand_marginal_prior(const blt_float_t* normal_lnprior,
                            const unsigned ref_gt,
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
    for (unsigned ngt(0); ngt<DIGT::SIZE; ++ngt)
    {
        grid_normal_lnprior[ngt] = (normal_lnprior[ngt]+ln_csse_rate);
    }

    // nostrand noise prior distributions for each allele combination axis:
    //
    // weight the prior by the potential originating genotypes:
    // if on AB axis, we want P(AA+noiseB)+P(AB+noise)+P(BB+noiseA)
    // so we have P(AA)*error_prob /3 + P(AB)*error_prob + P(BB)*error_prob/3
    //
    static const unsigned n_het_axes(6);
    blt_float_t nostrand_axis_prior[n_het_axes];
    for (unsigned ngt(N_BASE); ngt<DIGT::SIZE; ++ngt)
    {
        const unsigned axis_id(ngt-N_BASE);
        nostrand_axis_prior[axis_id] = normal_lnprior[ngt];
        // get the two associated homs:
        for (unsigned b(0); b<N_BASE; ++b)
        {
            if (DIGT::expect2(b,ngt)<=0) continue;
            nostrand_axis_prior[axis_id] = log_sum(nostrand_axis_prior[axis_id],
                                                   normal_lnprior[b]+ln_one_third);
        }
    }

    static const blt_float_t error_mod( -std::log(static_cast<blt_float_t>(DIGT_SGRID::HET_RES*2)) );

    // fill in normal sample prior for 'noise' frequencies:
    for (unsigned ngt(DIGT::SIZE); ngt<DIGT_SGRID::PRESTRAND_SIZE; ++ngt)
    {
        // 'ngt2' is the root diploid state corresponding to noise
        // state 'ngt'
        const unsigned ngt2(DIGT_SGRID::get_digt_state(ngt,ref_gt));
        assert(ngt2>=N_BASE);
        const unsigned axis_id(ngt2-N_BASE);
        grid_normal_lnprior[ngt] = (nostrand_axis_prior[axis_id]+ln_nostrand_sse_rate+error_mod);
        //        grid_normal_lnprior[ngt] = (normal_lnprior[ngt2]+ln_sse_rate+error_mod);
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
          const unsigned ref_gt,
          const blt_float_t sse_rate,
          const blt_float_t sseb_fraction,
          const blt_float_t somatic_normal_noise_rate,
          const bool is_somatic_normal_noise_rate,
          std::vector<blt_float_t>& grid_normal_lnprior,
          std::vector<blt_float_t>& somatic_marginal_lnprior)
{
    get_nostrand_marginal_prior(normal_lnprior,ref_gt,sse_rate,sseb_fraction,grid_normal_lnprior);
    if (is_somatic_normal_noise_rate)
    {
        get_nostrand_marginal_prior(normal_lnprior,ref_gt,somatic_normal_noise_rate,0,somatic_marginal_lnprior);
    }
    else
    {
        get_nostrand_marginal_prior(normal_lnprior,ref_gt,sse_rate,sseb_fraction,somatic_marginal_lnprior);
    }

    const blt_float_t strand_sse_rate(sse_rate*sseb_fraction);
    //    const blt_float_t nostrand_sse_rate(sse_rate-strand_sse_rate);

    //    const blt_float_t ln_csse_rate( log1p_switch(-sse_rate) );
    const blt_float_t ln_strand_sse_rate( std::log(strand_sse_rate) );
    //    const blt_float_t ln_nostrand_sse_rate( std::log(nostrand_sse_rate) );

    static const blt_float_t error_mod( -std::log(static_cast<blt_float_t>(DIGT_SGRID::HET_RES*2)) );

    // strand noise prior distributions for each allele combination axis:
    //
    // weight the prior by the potential originating genotypes: if on
    // AB axis, with A==ref, we want P(AA+noiseB) as the prior in the
    // model for the strand error states: the remaining term:
    // P(AB+noise)+P(BB+noiseA) is enumerated in to caluclate the
    // prior "throwaway state" -- regions of the stranded error
    // distribution for which we will approximate the lhood as 0
    //

    static const unsigned n_strand_het_axes(3);
    blt_float_t strand_axis_prior[n_strand_het_axes];

    for (unsigned sgt(0); sgt<n_strand_het_axes; ++sgt)
    {
        strand_axis_prior[sgt] = normal_lnprior[ref_gt]+ln_one_third;
    }

    // flaw in using error mod -- because strand state could exist and
    // be detectable at the canonical gt frequencies, we leave it the
    // same with the goal of stable performance as the user changes
    // the fraction term on the command-line:
    //
    // TODO: unclear comment?
    //

    for (unsigned ngt(DIGT_SGRID::PRESTRAND_SIZE); ngt<(DIGT_SGRID::SIZE); ++ngt)
    {
        const unsigned sgt(DIGT_SGRID::get_strand_state(ngt));
        grid_normal_lnprior[ngt] = (strand_axis_prior[sgt]+ln_strand_sse_rate+error_mod);
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
    _ln_som_mismatch=(std::log(opt.somatic_snv_rate/(static_cast<blt_float_t>((DIGT_SGRID::PRESTRAND_SIZE)-1))));

    for (unsigned i(0); i<(N_BASE+1); ++i)
    {
        prior_set& ps(_lnprior[i]);
        std::fill(ps.normal.begin(),ps.normal.end(),0);
        std::fill(ps.normal_poly.begin(),ps.normal_poly.end(),0);
    }

    for (unsigned i(0); i<(N_BASE+1); ++i)
    {
        prior_set& ps(_lnprior[i]);
        get_prior(pd_caller.lnprior_genomic(i),i,
                  opt.shared_site_error_rate,
                  opt.shared_site_error_strand_bias_fraction,
                  opt.site_somatic_normal_noise_rate,
                  opt.is_site_somatic_normal_noise_rate,
                  ps.normal,
                  ps.somatic_marginal);
        get_prior(pd_caller.lnprior_polymorphic(i),i,
                  opt.shared_site_error_rate,
                  opt.shared_site_error_strand_bias_fraction,
                  opt.site_somatic_normal_noise_rate,
                  opt.is_site_somatic_normal_noise_rate,
                  ps.normal_poly,
                  ps.somatic_marginal_poly);

        // special nostrand distro is used for somatic_gvcf:
        get_nostrand_marginal_prior(pd_caller.lnprior_genomic(i),i,opt.shared_site_error_rate,0,ps.normal_nostrand);
        get_nostrand_marginal_prior(pd_caller.lnprior_polymorphic(i),i,opt.shared_site_error_rate,0,ps.normal_poly_nostrand);
    }
}



// A simple static sized array with deep copy semantics:
//
template <unsigned NVAL>
struct cache_val
{
    std::array<blt_float_t,NVAL> val;
};



// This structure manages the caching of arrays of type:
// blt_float_t[NVAL] which are each associated with one value of
// (qscore,ratio_index).
//
// The get_val() function returns a tuple, the first value of which is
// a bool indicating whether the returned data structure has already
// been called, and thus is (presumably) cached. Note that the client is
// responsible for setting any values into the returned array for
// caching.
//
//
// this value caching didn't do much for the grid model -- better to
// leave it out for now... (how up to date is this comment?)
//
template <unsigned NVAL>
struct het_ratio_cache
{
    het_ratio_cache()
        : _is_cached(MAX_QSCORE* MAX_INDEX,false)
        , _cache(MAX_QSCORE* MAX_INDEX)
    {}

    std::pair<bool,cache_val<NVAL>*>
    get_val(const unsigned qscore,
            const unsigned ratio_index)
    {
        if (qscore>=MAX_QSCORE ||
            ratio_index>=MAX_INDEX)
        {
            return std::make_pair(false,&_any_val);
        }

        const unsigned index(ratio_index + qscore*MAX_INDEX);
        if (_is_cached[index])
        {
            return std::make_pair(true,&(_cache[index]));
        }
        else
        {
            _is_cached[index] = true;
            return std::make_pair(false,&(_cache[index]));
        }
    }

private:
    enum contanst { MAX_QSCORE = 64, MAX_INDEX = 12 };

    typedef cache_val<NVAL> cache_val_n;

    cache_val_n _any_val; // return this if a request is outside of cached range
    std::vector<bool> _is_cached;
    std::vector<cache_val_n> _cache;
};




// accelerated version with no hyrax q-val mods:
//
// the ratio key can be used as a proxy for the het ratio to look up cached results:
//
static
void
get_high_low_het_ratio_lhood_spi(const snp_pos_info& pi,
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
            const blt_float_t ceprob(1.-eprob);
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
increment_het_ratio_lhood_spi(const snp_pos_info& pi,
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
    get_high_low_het_ratio_lhood_spi(pi,het_ratio,het_ratio_index,hrcache,lhood_high,lhood_low);

    for (unsigned gt(0); gt<DIGT::SIZE; ++gt)
    {
        if (! DIGT::is_het(gt)) continue;
        all_het_lhood[gt] = log_sum(all_het_lhood[gt],lhood_high[gt]);
        all_het_lhood[gt] = log_sum(all_het_lhood[gt],lhood_low[gt]);
    }
}



// Fill in canonical dipliod positions in the likelihood
// function. This version is similar to the single sample version of
// the code, except that it doesn't allow any hyrax q-val adjustment
// shenanigans...
//
static
void
get_diploid_gt_lhood_spi(const blt_options& opt,
                         const snp_pos_info& pi,
                         const bool is_het_bias,
                         const blt_float_t het_bias,
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
            const blt_float_t ceprob(1.-eprob);
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
    if (is_het_bias)
    {
        // ! not thread-safe !
        static het_ratio_cache<3> hrcache_bias;

        // loop is currently setup to assume a uniform het ratio subgenotype prior
        const unsigned n_bias_steps(1+static_cast<unsigned>(het_bias/opt.het_bias_max_ratio_inc));
        const blt_float_t ratio_increment(het_bias/static_cast<blt_float_t>(n_bias_steps));
        for (unsigned i(0); i<n_bias_steps; ++i)
        {
            const blt_float_t het_ratio(0.5+(i+1)*ratio_increment);
            increment_het_ratio_lhood_spi(pi,het_ratio,i,hrcache_bias,lhood);
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
static
void
get_diploid_het_grid_lhood_spi(const snp_pos_info& pi,
                               blt_float_t* const lhood)
{
    // ! not thread-safe !
    static het_ratio_cache<3> hrcache;

    // get likelihood of each genotype
    for (unsigned gt(0); gt<(DIGT_SGRID::PRESTRAND_SIZE-DIGT::SIZE); ++gt) lhood[gt] = 0.;

    blt_float_t* lhood_off=lhood-N_BASE;

    static const blt_float_t ratio_increment(0.5/static_cast<blt_float_t>(DIGT_SGRID::HET_RES+1));
    for (unsigned i(0); i<DIGT_SGRID::HET_RES; ++i)
    {
        const blt_float_t het_ratio((i+1)*ratio_increment);
        get_high_low_het_ratio_lhood_spi(pi,het_ratio,i,hrcache,
                                         lhood_off+(i*DIGT_SGRID::HET_SIZE),
                                         lhood_off+((2*DIGT_SGRID::HET_RES-1-i)*DIGT_SGRID::HET_SIZE));
    }
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
                                   lhood+(i*DIGT_SGRID::STRAND_SIZE));
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
static
float
gvcf_nonsomatic_gvcf_prior(
    const somatic_snv_caller_strand_grid::prior_set& /*pset*/,
    const unsigned ngt,
    const unsigned tgt)
{
    if (ngt == tgt)
    {
        static const float lone(std::log(1.f));
        return lone;
    }
    else if (ngt<N_BASE)
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



// Given the likelihood, go through the final computations to get the
// posterior and derived values.
//
static
void
calculate_result_set_grid(
    const bool isComputeNonSomatic,
    const blt_float_t* normal_lhood,
    const blt_float_t* tumor_lhood,
    const somatic_snv_caller_strand_grid::prior_set& pset,
    const blt_float_t lnmatch,
    const blt_float_t lnmismatch,
    const unsigned /*ref_gt*/,
    const bool is_forced_output,
    result_set& rs)
{
    // a piece transplanted from 1150 to make a formal correction to
    // the priors which should have a low-impact on the results.
    // the prior below is incomplete
#ifdef DEBUG_ALTERNATE_PRIOR
    static const double neginf(-std::numeric_limits<double>::infinity());

    std::vector<double> prior(DDIGT_SGRID::SIZE);
    std::fill(prior.begin(),prior.end(),neginf);

    // this zero'd code is incomplete and abandoned for now...:
#if 0
    for (unsigned ngt(0); ngt<DIGT_SGRID::PRESTRAND_SIZE; ++ngt)
    {
        double base_prior(neginf);
        const bool is_noise(ngt>=STAR_DIINDEL::SIZE);
        if (is_noise)
        {
            base_prior=pset.normal[ngt];
        }
        else
        {
            base_prior=pset.nonoise[ngt];
        }
        for (unsigned tgt(0); tgt<DIGT_SGRID::PRESTRAND_SIZE; ++tgt)
        {
            const blt_float_t tgt_prior_mod( (tgt==ngt) ? lnmatch : lnmismatch );
            const unsigned dgt(DDIGT_SGRID::get_state(ngt,tgt));
            prior[dgt] = normal_genomic_lnprior[ngt]+tgt_prior_mod;
        }
    }

    for (unsigned gt(DIGT_SGRID::PRESTRAND_SIZE); gt<DIGT_SGRID::SIZE; ++gt)
    {
        const unsigned dgt(DDIGT_SGRID::get_state(gt,gt));
        prior[dgt] = normal_genomic_lnprior[gt]+lnmatch;
    }
#endif

    check_ln_distro(prior.begin(),
                    prior.end(),
                    "somatic snv full prior");
#endif

    // intentionally use higher float res (and heap alloc) for this structure:
    std::vector<double> pprob(DDIGT_SGRID::SIZE);

    // mult by prior distro to get unnormalized pprob for states in
    // the regular grid model:
    //
    for (unsigned ngt(0); ngt<DIGT_SGRID::PRESTRAND_SIZE; ++ngt)
    {
        for (unsigned tgt(0); tgt<DIGT_SGRID::PRESTRAND_SIZE; ++tgt)
        {
            const unsigned dgt(DDIGT_SGRID::get_state(ngt,tgt));

#if 0
            // the trusty old way...:
            const blt_float_t tgt_prior_mod( (tgt==ngt) ? lnmatch : lnmismatch );
            pprob[dgt] = normal_lhood[ngt]+tumor_lhood[tgt]+pset.normal[ngt]+tgt_prior_mod;
#else

            // unorm takes the role of the normal prior for the somatic case:
            //            static const blt_float_t unorm(std::log(static_cast<blt_float_t>(DIGT_SGRID::PRESTRAND_SIZE)));
            blt_float_t prior;
            if (tgt==ngt)
            {
                prior=pset.normal[ngt]+lnmatch;
            }
            else
            {
                prior=pset.somatic_marginal[ngt]+lnmismatch;
            }
            pprob[dgt] = normal_lhood[ngt]+tumor_lhood[tgt]+prior;
#endif
        }
    }

    // Now add the single-strand noise states. note that these states
    // are unique in that we don't look for mixtures of somatic
    // variation with these noise states, b/c single-strand
    // observations can almost exclusively be ruled out as noise:
    //
    for (unsigned gt(DIGT_SGRID::PRESTRAND_SIZE); gt<DIGT_SGRID::SIZE; ++gt)
    {
        const unsigned dgt(DDIGT_SGRID::get_state(gt,gt));
        pprob[dgt] = normal_lhood[gt]+tumor_lhood[gt]+pset.normal[gt]+lnmatch;
    }

    opt_normalize_ln_distro(pprob.begin(),pprob.end(),DDIGT_SGRID::is_nonsom.val.begin(),rs.max_gt);
    //normalize_ln_distro(pprob.begin(),pprob.end(),rs.max_gt);

    double nonsomatic_sum(0);
    for (unsigned gt(0); gt<DIGT_SGRID::SIZE; ++gt)
    {
        nonsomatic_sum += pprob[DDIGT_SGRID::get_state(gt,gt)];
    }
    rs.snv_qphred=error_prob_to_qphred(nonsomatic_sum);

    if ((! (is_forced_output || isComputeNonSomatic)) && (0==rs.snv_qphred)) return;

#if 0
    // alternate way to calculate the joint:
    //
    double min_not_somfrom_sum(0);
    for (unsigned dgt(0); dgt<DIGT::SIZE; ++dgt)
    {
        double not_somfrom_sum(nonsomatic_sum);

        for (unsigned ngt(0); ngt<DIGT_SGRID::PRESTRAND_SIZE; ++ngt)
        {
            // we're looking for the joint prob when state dgt is true
            // in the normal, so skip this as a normal state here:
            //
            if (dgt==ngt) continue;

            for (unsigned tgt(0); tgt<DIGT_SGRID::PRESTRAND_SIZE; ++tgt)
            {
                // we've already started from the nonsomatic som, so we can skip the equal states:
                //
                if (ngt==tgt) continue;

                not_somfrom_sum += pprob[DDIGT_SGRID::get_state(ngt,tgt)];
            }
        }

        if ((dgt==0) || (!_somfrom_sum<min_not_somfrom_sum))
        {
            min_not_somfrom_sum=not_somfrom_sum;
            rs.snv_from_ntype_qphred=error_prob_to_qphred(not_somfrom_sum);
            rs.ntype=dgt;
        }
    }
#endif

#if 0
    // reset max_gt to the most likely state excluding normal noise states:
    //
    rs.max_gt=0;
    for (unsigned dgt(0); dgt<DIGT::SIZE; ++dgt)
    {
        for (unsigned tgt(0); tgt<DIGT_SGRID::PRESTRAND_SIZE; ++tgt)
        {
            const unsigned xgt(DDIGT_SGRID::get_state(dgt,tgt));
            if (pprob[xgt] > pprob[rs.max_gt]) rs.max_gt=xgt;
        }
    }
#endif

    // Calculate normal distribution alone so that we can classify this call:
    //
    // Polymorphic prior is used because in this situation we want to
    // be conservative about the reference classification --
    // ie. conditioned on only looking at putative somatic sites, we
    // require evidence to show that the normal is in fact reference
    // and not simply an unsampled copy of the somatic variation.
    //
    std::vector<double> normal_pprob(DIGT_SGRID::PRESTRAND_SIZE);
    for (unsigned ngt(0); ngt<DIGT_SGRID::PRESTRAND_SIZE; ++ngt)
    {
        normal_pprob[ngt] = normal_lhood[ngt]+pset.normal_poly[ngt];
    }

    unsigned max_norm_gt(0);
    normalize_ln_distro(normal_pprob.begin(),normal_pprob.end(),max_norm_gt);

    // find the probability of max_norm_gt:
    const double ngt_prob(prob_comp(normal_pprob.begin(),normal_pprob.end(),max_norm_gt));

    // (1-(1-a)(1-b)) -> a+b-(ab)
    double not_somfrom_sum(nonsomatic_sum+ngt_prob-(nonsomatic_sum*ngt_prob));

    rs.snv_from_ntype_qphred=error_prob_to_qphred(not_somfrom_sum);
    rs.ntype=max_norm_gt;

    // add new somatic gVCF value -- note this is an expanded definition of 'non-somatic' beyond just f_N == f_T
    if (isComputeNonSomatic)
    {
        /// process regular tumor/normal lhood, but:
        /// (1) use uniform probability for {somatic,non-somatic} states
        /// (2) simply computation to remove strand-specific logic
        /// (3) ignore normal genotype
        ///
        for (unsigned ngt(0); ngt<DIGT_SGRID::PRESTRAND_SIZE; ++ngt)
        {
            for (unsigned tgt(0); tgt<DIGT_SGRID::PRESTRAND_SIZE; ++tgt)
            {
                const unsigned dgt(DDIGT_SGRID::get_state(ngt,tgt));
                pprob[dgt] = normal_lhood[ngt]+tumor_lhood[tgt]+gvcf_nonsomatic_gvcf_prior(pset,ngt,tgt);
            }
        }

        unsigned max_gt(0);
        opt_normalize_ln_distro(pprob.begin(),pprob.begin()+DDIGT_SGRID::PRESTRAND_SIZE,
                                DDIGT_SGRID::is_nonsom.val.begin(),max_gt);

        double sgvcf_nonsomatic_sum(0);
        for (unsigned ngt(0); ngt<DIGT_SGRID::PRESTRAND_SIZE; ++ngt)
        {
            for (unsigned tgt(0); tgt<DIGT_SGRID::PRESTRAND_SIZE; ++tgt)
            {
                if (ngt != tgt) continue;
                const unsigned dgt(DDIGT_SGRID::get_state(ngt,tgt));
                sgvcf_nonsomatic_sum += pprob[dgt];
            }
        }

        rs.nonsomatic_qphred=error_prob_to_qphred(1.-sgvcf_nonsomatic_sum);
    }

    static const bool is_compute_sb(true);
    if (is_compute_sb)
    {
        // get ratio of strand bias vs. non-strand-bias version of max_gt, if max_gt does not correspond to a het state, then
        // set sb to 0
        unsigned normal_gt,tumor_gt;
        DDIGT_SGRID::get_digt_grid_states(
            rs.max_gt,
            normal_gt,
            tumor_gt);

//        const unsigned het_count(DIGT_SGRID::get_het_count(tumor_gt));

//        if ((tumor_gt>=N_BASE) && (het_count < DIGT_SGRID::STRAND_COUNT))
        {
#if 0
            const bool is_strand_state(DIGT_SGRID::is_strand_state(tumor_gt));

            unsigned symm_tumor_gt(0);
            unsigned strand_tumor_gt(0);
            if (! is_strand_state)
            {
                strand_tumor_gt=DIGT_SGRID::toggle_strand_state(tumor_gt, ref_base);
                symm_tumor_gt=tumor_gt;
            }
            else
            {
                strand_tumor_gt=tumor_gt;
                symm_tumor_gt=DIGT_SGRID::toggle_strand_state(tumor_gt, ref_base);
            }
            rs.strandBias = tumor_lhood[strand_tumor_gt] - tumor_lhood[symm_tumor_gt];
#endif

            const blt_float_t symm_lhood(*std::max_element(tumor_lhood+N_BASE, tumor_lhood+DIGT_SGRID::PRESTRAND_SIZE));
            const blt_float_t strand_lhood(*std::max_element(tumor_lhood+DIGT_SGRID::PRESTRAND_SIZE, tumor_lhood+DIGT_SGRID::SIZE));
            rs.strandBias = std::max(0.f,(strand_lhood - symm_lhood));
        }
#if 0
        else
        {
            rs.strandBias = 0.;
        }
#endif
    }
}



///
///
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

        if (normal_pi.ref_base=='N')
        {
            sgt.is_forced_output=false;
            return;
        }
        sgt.ref_gt=base_to_id(normal_pi.ref_base);

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

        // get likelihood of each genotype
        //
        static constexpr bool is_normal_het_bias(false);
        static constexpr blt_float_t normal_het_bias(0.0);
        static constexpr bool is_tumor_het_bias(false);
        static constexpr blt_float_t tumor_het_bias(0.0);

        const extended_pos_info& nepi(is_include_tier2 ? *normal_epi_t2_ptr : normal_epi );
        const extended_pos_info& tepi(is_include_tier2 ? *tumor_epi_t2_ptr : tumor_epi );
        get_diploid_gt_lhood_spi(_opt,nepi.pi,is_normal_het_bias,normal_het_bias,normal_lhood);
        get_diploid_gt_lhood_spi(_opt,tepi.pi,is_tumor_het_bias,tumor_het_bias,tumor_lhood);

        get_diploid_het_grid_lhood_spi(nepi.pi,normal_lhood+DIGT::SIZE);
        get_diploid_het_grid_lhood_spi(tepi.pi,tumor_lhood+DIGT::SIZE);

        get_diploid_strand_grid_lhood_spi(nepi.pi,sgt.ref_gt,normal_lhood+DIGT_SGRID::PRESTRAND_SIZE);
        get_diploid_strand_grid_lhood_spi(tepi.pi,sgt.ref_gt,tumor_lhood+DIGT_SGRID::PRESTRAND_SIZE);

        // genomic site results:
        calculate_result_set_grid(isComputeNonSomatic,
                                  normal_lhood,
                                  tumor_lhood,
                                  get_prior_set(sgt.ref_gt),
                                  _ln_som_match,_ln_som_mismatch,
                                  sgt.ref_gt,
                                  sgt.is_forced_output,
                                  tier_rs[i]);

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
        if       (sgt.rs.ntype==sgt.ref_gt)
        {
            sgt.rs.ntype=NTYPE::REF;
        }
        else if (DIGT::is_het(sgt.rs.ntype))
        {
            sgt.rs.ntype=NTYPE::HET;
        }
        else
        {
            sgt.rs.ntype=NTYPE::HOM;
        }
    }

    sgt.rs.snv_qphred = tier_rs[sgt.snv_tier].snv_qphred;

    /// somatic gVCF, always use tier1 to keep things simple:
    sgt.rs.nonsomatic_qphred = tier_rs[0].nonsomatic_qphred;
}
