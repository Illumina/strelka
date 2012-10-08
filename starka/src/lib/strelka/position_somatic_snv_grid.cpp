// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
//

/// \file

/// \author Chris Saunders
///

#if 0

#include "position_somatic_snv_grid.hh"

#include "blt_common/snp_util.hh"
#include "blt_util/log.hh"
#include "blt_util/math_util.hh"
#include "blt_util/prob_util.hh"
#include "blt_util/seq_util.hh"

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <iostream>
#include <iomanip>
#include <map>


const blt_float_t one_third(1./3.);
const blt_float_t ln_one_third(std::log(one_third));
const blt_float_t one_half(1./2.);
const blt_float_t ln_one_half(std::log(one_half));



namespace DDIGT_GRID {

    is_nonsom_maker_t::
    is_nonsom_maker_t()
      : val(SIZE,false)
    {
        for(unsigned gt(0);gt<DIGT_GRID::SIZE;++gt){
            val[get_state(gt,gt)] = true;
        }
    }   

    const is_nonsom_maker_t is_nonsom;
}





std::ostream&
operator<<(std::ostream& os,
           const DDIGT_GRID::index_t dgt) {

    unsigned normal_gt;
    unsigned tumor_gt;
    DDIGT_GRID::get_digt_grid_states(dgt,normal_gt,tumor_gt);

    os << DIGT::label(DIGT_GRID::get_digt_state(normal_gt)) 
       << "->"
       << DIGT::label(DIGT_GRID::get_digt_state(tumor_gt));

    return os;

}


#if 0
static
void
get_prior(const blt_float_t* normal_lnprior,
          std::vector<blt_float_t>& grid_normal_lnprior) {
    
    static const blt_float_t het_mod( -std::log(static_cast<blt_float_t>(DIGT_GRID::HET_COUNT)) );
    
    for(unsigned ngt(0);ngt<DIGT_GRID::SIZE;++ngt) {
        const unsigned ngt2(DIGT_GRID::get_digt_state(ngt));
        const bool is_het(DIGT::is_het(ngt2));
        const blt_float_t normal_mod(is_het ? het_mod : 0 );
        grid_normal_lnprior[ngt] = (normal_lnprior[ngt2]+normal_mod);
    }
}
#endif



// ln_sse_rate = log( P(S) )
// ln_csse_rate = log( 1-P(S) )
//
// grid_normal_lnprior is the mixture of normal diploid probabilities with uniform noise
//
static
void
get_prior(const blt_float_t* normal_lnprior,
          const blt_float_t sse_rate,
          std::vector<blt_float_t>& grid_normal_lnprior) {

    static const blt_float_t error_mod( -std::log(static_cast<blt_float_t>(DIGT_GRID::HET_RES*2)) );

    const blt_float_t ln_sse_rate( std::log(sse_rate) );
    const blt_float_t ln_csse_rate( log1p_switch(-sse_rate) );

    for(unsigned ngt(0);ngt<DIGT::SIZE;++ngt) {
        grid_normal_lnprior[ngt] = (normal_lnprior[ngt]+ln_csse_rate);
    }

    // weight the prior by the potential originating genotypes:
    // if on AB axis, we want P(AA+noiseB)+P(AB+noise)+P(BB+noiseA)
    // so we have P(AA)*error_prob /3 + P(AB)*error_prob + P(BB)*error_prob/3
    //
    static const unsigned n_het_axes(6);
    blt_float_t axis_prior[n_het_axes];
    for(unsigned ngt(N_BASE);ngt<DIGT::SIZE;++ngt) {
        const unsigned axis_id(ngt-N_BASE);
        axis_prior[axis_id] = normal_lnprior[ngt];
        // get the two associated homs:
        for(unsigned b(0);b<N_BASE;++b){
            if(DIGT::expect2(b,ngt)<=0) continue;
            axis_prior[axis_id] = log_sum(axis_prior[axis_id],
                                          normal_lnprior[b]+ln_one_third);
        }
    }

    for(unsigned ngt(DIGT::SIZE);ngt<DIGT_GRID::SIZE;++ngt) {
        // ngt2 indicates the "root" heterozygote state of the grid state:
        const unsigned ngt2(DIGT_GRID::get_digt_state(ngt));
        assert(ngt2>=N_BASE);
        const unsigned axis_id(ngt2-N_BASE);
        grid_normal_lnprior[ngt] = (axis_prior[axis_id]+ln_sse_rate+error_mod);
        //        grid_normal_lnprior[ngt] = (normal_lnprior[ngt2]+ln_sse_rate+error_mod);
    }

#ifdef SOMATIC_DEBUG
    check_ln_distro(grid_normal_lnprior.begin(),
                    grid_normal_lnprior.end(),
                    "normal somatic prior");
#endif
}



somatic_snv_caller_grid::
somatic_snv_caller_grid(const strelka_options& opt,
                        const pprob_digt_caller& pd_caller) 
    : _opt(opt), _pd_caller(pd_caller) {

    _ln_som_match=(log1p_switch(-opt.somatic_snv_rate));
    _ln_som_mismatch=(std::log(opt.somatic_snv_rate/(static_cast<blt_float_t>((DIGT_GRID::SIZE)-1))));

    for(unsigned i(0);i<(N_BASE+1);++i) {
        prior_set& ps(_lnprior[i]);
        std::fill(ps.normal.begin(),ps.normal.end(),0);
        std::fill(ps.normal_poly.begin(),ps.normal_poly.end(),0);
    }

    for(unsigned i(0);i<(N_BASE+1);++i) {
        prior_set& ps(_lnprior[i]);
        get_prior(pd_caller.lnprior_genomic(i),
                  opt.shared_site_error_rate,
                  ps.normal);
        get_prior(pd_caller.lnprior_polymorphic(i),
                  opt.shared_site_error_rate,
                  ps.normal_poly);
    }

#ifdef SOMATIC_DEBUG
    std::vector<blt_float_t>& grid_normal_lnprior(_lnprior[0].normal);

    // double-check full prior:
    std::vector<blt_float_t> prior(DDIGT_GRID::SIZE);

    for(unsigned ngt(0);ngt<DIGT_GRID::SIZE;++ngt){
        for(unsigned tgt(0);tgt<DIGT_GRID::SIZE;++tgt){
            const blt_float_t tgt_prior_mod( (tgt==ngt) ? _ln_som_match : _ln_som_mismatch );
            const unsigned dgt(DDIGT_GRID::get_state(ngt,tgt));
            prior[dgt] = grid_normal_lnprior[ngt]+tgt_prior_mod;
        }
    }

    check_ln_distro(prior.begin(),
                    prior.end(),
                    "full somatic prior");
#endif
}


struct cache_val {

    enum { SIZE = 3 };

    cache_val() {}

    cache_val(const cache_val& x) {
        for(unsigned i(0);i<SIZE;++i){
            val[i] = x.val[i];
        }
    }

    cache_val&
    operator==(const cache_val& rhs){
        if(&rhs==this) return *this;
        for(unsigned i(0);i<SIZE;++i){
            val[i] = rhs.val[i];
        }
        return *this;
    }

    blt_float_t
    high_val(const unsigned i) const {
        return val[i];
    }

    blt_float_t
    low_val(const unsigned i) const {
        static const uint8_t remap[SIZE] = {0,2,1};
        return val[remap[i]];
    }

    blt_float_t val[SIZE];
};



// this value caching didn't do much for the grid model -- better to leave it out for now...
//
#if 0
struct het_ratio_cache {

    std::pair<bool,cache_val*>
    get_val(const unsigned qscore,
            const unsigned ratio_index) {

        map_t::iterator i(_rcache.find(std::make_pair(qscore,ratio_index)));
        if(i != _rcache.end()) return std::make_pair(true,&(i->second));
        return std::make_pair(false,&(_rcache[std::make_pair(qscore,ratio_index)]));
    }

private:
    typedef std::pair<unsigned,unsigned> key_t; // qscore,ratio_key
    typedef std::map<key_t,cache_val> map_t;

    map_t _rcache;
};



// no thread safety needed for now...
//
het_ratio_cache hrcache;
#endif



// accelerated version with no hyrax q-val mods:
//
// the ratio key can be used as a proxy for the het ratio to look up cached results:
//
static
void
get_high_low_het_ratio_lhood_spi(const snp_pos_info& pi,
                                 const blt_float_t het_ratio,
                                 const unsigned het_ratio_index,
                                 blt_float_t* lhood_high,
                                 blt_float_t* lhood_low) {

    const blt_float_t chet_ratio(1.-het_ratio);

    const unsigned n_calls(pi.calls.size());

    cache_val cv;

    for(unsigned i(0);i<n_calls;++i){
        const base_call& bc(pi.calls[i]);

#if 0
        std::pair<bool,cache_val*> ret(hrcache.get_val(bc.qscore,het_ratio_index));
        cache_val& cv(*ret.second);
        if(not ret.first)   
#endif
        const blt_float_t eprob(bc.error_prob());
        const blt_float_t ceprob(1.-eprob);
        //const blt_float_t lne(bc.ln_error_prob());
        //const blt_float_t lnce(bc.ln_comp_error_prob());

        // precalculate the result for expect values of 0.0, het_ratio, chet_ratio, 1.0
        cv.val[0] = bc.ln_error_prob()+ln_one_third;
        cv.val[1] = std::log((ceprob)*het_ratio+((eprob)*one_third)*chet_ratio);
        cv.val[2] = std::log((ceprob)*chet_ratio+((eprob)*one_third)*het_ratio);
        
        const uint8_t obs_id(bc.base_id);

        for(unsigned gt(N_BASE);gt<DIGT::SIZE;++gt){
            const unsigned key(DIGT::expect2_bias(obs_id,gt));
            lhood_high[gt] += cv.high_val(key);
            lhood_low[gt] += cv.low_val(key);
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
                              blt_float_t* all_het_lhood) {

    // multiply probs of alternate ratios into local likelihoods, then
    // *add* them to the global tally (effictively this is the sum lhood of
    // many different heterozygous genotypes).
    //
    // in the gt_high genotype, the first allele (in lexicographical
    // order) is expected at het_ratio and the second allele is
    // expected at chet_ratio.  gt_low genotype is vice versa.
    //
    blt_float_t lhood_high[DIGT::SIZE];
    blt_float_t lhood_low[DIGT::SIZE];
    for(unsigned gt(0);gt<DIGT::SIZE;++gt) {
        lhood_high[gt] = 0.;
        lhood_low[gt] = 0.;
    }
    get_high_low_het_ratio_lhood_spi(pi,het_ratio,het_ratio_index,lhood_high,lhood_low);

    for(unsigned gt(0);gt<DIGT::SIZE;++gt){
        if(not DIGT::is_het(gt)) continue;
        all_het_lhood[gt] = log_sum(all_het_lhood[gt],lhood_high[gt]);
        all_het_lhood[gt] = log_sum(all_het_lhood[gt],lhood_low[gt]);
    }
}



// this version doesn't allow any hyrax q-val shenanigans...
//
static
void
get_diploid_gt_lhood_spi(const blt_options& opt,
                         const snp_pos_info& pi,
                         const bool is_het_bias,
                         const blt_float_t het_bias,
                         blt_float_t* const lhood) {

    // get likelihood of each genotype
    for(unsigned gt(0);gt<DIGT::SIZE;++gt) lhood[gt] = 0.;

    const unsigned n_calls(pi.calls.size());
    for(unsigned i(0);i<n_calls;++i){
        const base_call& bc(pi.calls[i]);
        const blt_float_t eprob(bc.error_prob());
        const blt_float_t ceprob(1.-eprob);
        const blt_float_t lne(bc.ln_error_prob());
        const blt_float_t lnce(bc.ln_comp_error_prob());

        // precalculate the result for expect values of 0.0, 0.5 & 1.0
        blt_float_t val[3];
        val[0] = lne+ln_one_third;
        val[1] = std::log((ceprob)+((eprob)*one_third))+ln_one_half;
        val[2] = lnce;

        const uint8_t obs_id(pi.calls[i].base_id);
        for(unsigned gt(0);gt<DIGT::SIZE;++gt){
            lhood[gt] += val[DIGT::expect2(obs_id,gt)];
        }
    }

    if(is_het_bias) {
        // loop is currently setup to assume a uniform het ratio subgenotype prior
        const unsigned n_bias_steps(1+static_cast<unsigned>(het_bias/opt.het_bias_max_ratio_inc));
        const blt_float_t ratio_increment(het_bias/static_cast<blt_float_t>(n_bias_steps));
        for(unsigned i(0);i<n_bias_steps;++i) {
            const blt_float_t het_ratio(0.5+(i+1)*ratio_increment);
            increment_het_ratio_lhood_spi(pi,het_ratio,i,lhood);
        }
        
        const unsigned n_het_subgt(1+2*n_bias_steps);
        const blt_float_t subgt_log_prior(std::log(1./static_cast<blt_float_t>(n_het_subgt)));
        
        for(unsigned gt(0);gt<DIGT::SIZE;++gt){
            if(not DIGT::is_het(gt)) continue;
            lhood[gt] += subgt_log_prior;
        }
    }


}



static
void
get_diploid_het_grid_lhood_spi(const snp_pos_info& pi,
                               blt_float_t* const lhood) {

    // get likelihood of each genotype
    for(unsigned gt(0);gt<(DIGT_GRID::SIZE-DIGT::SIZE);++gt) lhood[gt] = 0.;

    blt_float_t* lhood_off=lhood-N_BASE;

    static const blt_float_t ratio_increment(0.5/static_cast<blt_float_t>(DIGT_GRID::HET_RES+1));
    for(unsigned i(0);i<DIGT_GRID::HET_RES;++i) {
        const blt_float_t het_ratio((i+1)*ratio_increment);
        get_high_low_het_ratio_lhood_spi(pi,het_ratio,i,
                                         lhood_off+(i*DIGT_GRID::HET_SIZE),
                                         lhood_off+((2*DIGT_GRID::HET_RES-1-i)*DIGT_GRID::HET_SIZE));
    }
}



typedef somatic_snv_genotype_grid::result_set result_set;


#ifdef SOMATIC_DEBUG
static
void
debug_dump_ddigt_lhood(const blt_float_t* lhood,
                       std::ostream& os) {

    blt_float_t pprob[DDIGT::SIZE];
    for(unsigned gt(0);gt<DDIGT::SIZE;++gt){
        pprob[gt] = lhood[gt];
    }

    unsigned max_gt(0);
    normalize_ln_prob(pprob,pprob+DDIGT::SIZE,max_gt);

    os << std::setprecision(3) << std::fixed;
    for(unsigned ngt(0);ngt<DIGT::SIZE;++ngt){
        for(unsigned tgt(0);tgt<DIGT::SIZE;++tgt){
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
calculate_result_set_grid(const blt_float_t* normal_lhood,
                          const blt_float_t* tumor_lhood,
                          const std::vector<blt_float_t>& normal_genomic_lnprior,
                          const std::vector<blt_float_t>& normal_polymorphic_lnprior,
                          const blt_float_t lnmatch,
                          const blt_float_t lnmismatch,
                          const unsigned ref_gt,
                          result_set& rs) {

    std::vector<blt_float_t> pprob(DDIGT_GRID::SIZE);

    // mult by prior distro to get unnormalized pprob:
    //
    for(unsigned ngt(0);ngt<DIGT_GRID::SIZE;++ngt){
        //        const bool is_ngt_het(DIGT::is_het(DIGT_GRID::get_digt_state(ngt)));
        for(unsigned tgt(0);tgt<DIGT_GRID::SIZE;++tgt){
            //  const bool is_tgt_het(DIGT::is_het(tgt));
            const blt_float_t tgt_prior_mod( (tgt==ngt) ? lnmatch : lnmismatch );
            const unsigned dgt(DDIGT_GRID::get_state(ngt,tgt));
            pprob[dgt] = normal_lhood[ngt]+tumor_lhood[tgt]+normal_genomic_lnprior[ngt]+tgt_prior_mod;
        }
    }

    opt_normalize_ln_distro(pprob.begin(),pprob.end(),DDIGT_GRID::is_nonsom.val.begin(),rs.max_gt);
    //normalize_ln_prob(pprob.begin(),pprob.end(),rs.max_gt);

    blt_float_t nonsomatic_sum(0);
    for(unsigned gt(0);gt<DIGT_GRID::SIZE;++gt){
        nonsomatic_sum += pprob[DDIGT_GRID::get_state(gt,gt)];
    }
    rs.snv_qphred=error_prob_to_qphred(nonsomatic_sum);

    if(0==rs.snv_qphred) return;

    // calculate normal distribution alone so that we can classify this call:
    std::vector<blt_float_t> normal_pprob(DIGT_GRID::SIZE);
    for(unsigned ngt(0);ngt<DIGT_GRID::SIZE;++ngt){
        normal_pprob[ngt] = normal_lhood[ngt]+normal_polymorphic_lnprior[ngt];
    }

    unsigned max_norm_gt(0);
    normalize_ln_distro(normal_pprob.begin(),normal_pprob.end(),max_norm_gt);

    const blt_float_t nonref_prob(prob_comp(normal_pprob.begin(),normal_pprob.end(),ref_gt));

    // (1-(1-a)(1-b)) -> a+b-(ab)
    blt_float_t not_somfrom_ref_sum(nonsomatic_sum+nonref_prob-(nonsomatic_sum*nonref_prob));

    rs.snv_from_ref_qphred=error_prob_to_qphred(not_somfrom_ref_sum);
    
    rs.max_gt_qphred=error_prob_to_qphred(prob_comp(pprob.begin(),pprob.end(),rs.max_gt));
}



///
///
void
somatic_snv_caller_grid::
position_somatic_snv_call(const extended_pos_info& normal_epi,
                          const extended_pos_info& tumor_epi,
                          const extended_pos_info* normal_epi_t2_ptr,
                          const extended_pos_info* tumor_epi_t2_ptr,
                          somatic_snv_genotype_grid& sgt) const {

    static const bool is_always_test(false);

    {
        const snp_pos_info& normal_pi(normal_epi.pi);
        const snp_pos_info& tumor_pi(tumor_epi.pi);

        if(normal_pi.ref_base=='N') return;
        sgt.ref_gt=base_to_id(normal_pi.ref_base);

        // check that a non-reference call meeting quality criteria even
        // exists:
        if(not is_always_test) {
            if(is_spi_allref(normal_pi,sgt.ref_gt) and is_spi_allref(tumor_pi,sgt.ref_gt)) return;
        }
    }

    // strawman model treats normal and tumor as independent, so
    // calculate separate lhoods:
    blt_float_t normal_lhood[DIGT_GRID::SIZE];
    blt_float_t tumor_lhood[DIGT_GRID::SIZE];

    const bool is_tier2(NULL != normal_epi_t2_ptr);

    static const unsigned n_tier(2);
    result_set tier_rs[n_tier];
    for(unsigned i(0);i<n_tier;++i) {
        const bool is_include_tier2(i==1);
        if(is_include_tier2){
            if(not is_tier2) continue;
            if(tier_rs[0].snv_qphred==0){
                tier_rs[1].snv_qphred=0;
                continue;
            }
        }

        // get likelihood of each genotype
        //
        static const bool is_normal_het_bias(false);
        static const blt_float_t normal_het_bias(0.0);
        static const bool is_tumor_het_bias(false);
        static const blt_float_t tumor_het_bias(0.0);

        const extended_pos_info& nepi(is_include_tier2 ? *normal_epi_t2_ptr : normal_epi );
        const extended_pos_info& tepi(is_include_tier2 ? *tumor_epi_t2_ptr : tumor_epi );
        get_diploid_gt_lhood_spi(_opt,nepi.pi,is_normal_het_bias,normal_het_bias,normal_lhood);
        get_diploid_gt_lhood_spi(_opt,tepi.pi,is_tumor_het_bias,tumor_het_bias,tumor_lhood);

        get_diploid_het_grid_lhood_spi(nepi.pi,normal_lhood+DIGT::SIZE);
        get_diploid_het_grid_lhood_spi(tepi.pi,tumor_lhood+DIGT::SIZE);

#if 0
        // get normal diploid results (will use these for somatic call classification):
        diploid_genotype::result_set normal_digt_rs;
        pprob_digt_caller::calculate_result_set(normal_lhood,
                                                _pd_caller.lnprior_genomic(sgt.ref_gt),
                                                sgt.ref_gt,
                                                normal_digt_rs);
#endif

        // genomic site results:
        calculate_result_set_grid(normal_lhood,
                                  tumor_lhood,
                                  lnprior_genomic(sgt.ref_gt),
                                  lnprior_polymorphic(sgt.ref_gt),
                                  _ln_som_match,_ln_som_mismatch,
                                  sgt.ref_gt,
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
    }

    sgt.tier=0;
    if(is_tier2) {

#if 1
    if(tier_rs[0].snv_from_ref_qphred > tier_rs[1].snv_from_ref_qphred){
            sgt.tier=1;
        }
#else
        if(tier_rs[0].snv_qphred > tier_rs[1].snv_qphred) {
            sgt.tier=1;
        }
#endif

    }
    sgt.genome=tier_rs[sgt.tier];


    // not sure which prior to use yet:
    //
    sgt.is_snv=((sgt.genome.snv_qphred != 0));
    //    sgt.is_snv=((sgt.genome.snv_qphred != 0) or (sgt.poly.snv_qphred != 0));


#ifdef SOMATIC_DEBUG
    if(sgt.is_snv) {
    log_os << "tier1_qphred: " << tier_rs[0].snv_qphred << "\n";
    log_os << "tier2_qphred: " << tier_rs[1].snv_qphred << "\n";
    result_set rs(sgt.genome);
    const blt_float_t* lnprior(lnprior_genomic(sgt.ref_gt));
    for(unsigned ngt(0);ngt<DIGT::SIZE;++ngt){
        for(unsigned tgt(0);tgt<DIGT::SIZE;++tgt){
            const unsigned dgt(DDIGT::get_state(ngt,tgt));
            rs.pprob[dgt] = normal_lhood[ngt]+tumor_lhood[tgt]+lnprior[dgt];
        }
    }

    log_os << "normal_digt: ";
    debug_dump_digt_lhood(normal_lhood,log_os);
    log_os << "\n";
    log_os << "tumor_digt: ";
    debug_dump_digt_lhood(tumor_lhood,log_os);
    log_os << "\n";
    log_os << "prior:\n";
    debug_dump_ddigt_lhood(lnprior,log_os);
    log_os << "pprob:\n";
    blt_float_t pprob[DDIGT::SIZE];
    for(unsigned gt(0);gt<DDIGT::SIZE;++gt){
        pprob[gt] = rs.pprob[gt];
    }   
    debug_dump_ddigt_lhood(pprob,log_os);
    }
#endif

}



static
void
write_result_set(const result_set& rs,
                 std::ostream& os) {

    os << rs.snv_qphred
       << '\t' << rs.snv_from_ref_qphred
       << '\t' << rs.snv_from_het_qphred
       << '\t' << rs.snv_from_het_loh_qphred
       << '\t' << rs.snv_from_het_nonloh_qphred
       << '\t' << rs.snv_from_hom_qphred
       << '\t' << rs.snv_from_anyhom_qphred
       << '\t' << static_cast<DDIGT_GRID::index_t>(rs.max_gt) 
       << '\t' << rs.max_gt_qphred;
}



void
write_somatic_snv_genotype_grid(const strelka_options& opt,
                                const somatic_snv_genotype_grid& sgt,
                                const snp_pos_info& normal_pi,
                                const snp_pos_info& tumor_pi,
                                std::ostream& os) {

    os << std::setprecision(10) << std::fixed;
    
    const result_set& ge(sgt.genome);
#ifdef ENABLE_POLY
    const result_set& po(sgt.poly);
#endif

    os << (sgt.tier+1) << '\t';
    write_result_set(ge,os);
#ifdef ENALBLE_POLY
    os << '\t';
    write_result_set(po,os);
#endif
    if(opt.is_print_used_allele_counts) {
        normal_pi.print_known_counts(os);
        tumor_pi.print_known_counts(os);
    }

     os.unsetf(std::ios::fixed);
}

#endif
