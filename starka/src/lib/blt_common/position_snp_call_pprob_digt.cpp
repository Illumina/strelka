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
#include "blt_common/position_snp_call_pprob_digt.hh"

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



const blt_float_t one_third(1./3.);
const blt_float_t log_one_third(std::log(one_third));
const blt_float_t one_half(1./2.);
const blt_float_t log_one_half(std::log(one_half));



static
void
get_genomic_prior(const unsigned ref_gt,
                  const blt_float_t theta,
                  blt_float_t* const prior) {

    blt_float_t prior_sum(0.);
    for(unsigned gt(0);gt<DIGT::SIZE;++gt){
        if(gt==ref_gt) continue;
        prior[gt]=(theta*one_third);
        if(DIGT::is_het(gt)){
            if(DIGT::expect(ref_gt,gt)<=0.) prior[gt]*=theta;
        } else {
            prior[gt]*=.5;
        }
        prior_sum += prior[gt];
    }
    assert(prior_sum <= 1.);
    prior[ref_gt] = (1.-prior_sum);
}



static
void
get_poly_prior(const unsigned ref_gt,
               const blt_float_t theta,
               blt_float_t* const prior) {

    blt_float_t prior_sum(0.);
    const blt_float_t ctheta(1.-theta);
    for(unsigned gt(0);gt<DIGT::SIZE;++gt){
        if(gt==ref_gt) {
            prior[gt]=0.25*(ctheta);
        } else if(DIGT::is_het(gt)){
            if(DIGT::expect(ref_gt,gt)<=0.) {
                prior[gt] = theta*one_third;
            } else {
                prior[gt] = 0.5*one_third*ctheta;
            }
        } else {
            prior[gt] = 0.25*one_third*ctheta;
        }
        prior_sum += prior[gt];
    }
    assert(std::abs(1.-prior_sum) < 0.0001);
}



static
void
sum_gt(blt_float_t* const x1,
       const blt_float_t* const x2) {
    for(unsigned gt(0);gt<DIGT::SIZE;++gt){ x1[gt] += x2[gt]; }
}



static
void
norm_gt(blt_float_t* const x) {
    blt_float_t sum(0);
    for(unsigned gt(0);gt<DIGT::SIZE;++gt){ sum += x[gt]; }
    sum = 1./sum;
    for(unsigned gt(0);gt<DIGT::SIZE;++gt){ x[gt] *= sum; }
}



static
void
log_gt(blt_float_t* const x) {
    for(unsigned gt(0);gt<DIGT::SIZE;++gt){ x[gt] = std::log(x[gt]); }
}



pprob_digt_caller::
pprob_digt_caller(const blt_float_t theta) {

    for(unsigned i(0);i<(N_BASE+1);++i) {
        prior_set& ps(_lnprior[i]);
        std::fill(ps.genome,ps.genome+DIGT::SIZE,0);
        std::fill(ps.poly,ps.poly+DIGT::SIZE,0);
    }

    for(unsigned i(0);i<N_BASE;++i) {
        prior_set& ps(_lnprior[i]);
        get_genomic_prior(i,theta,ps.genome);
        get_poly_prior(i,theta,ps.poly);
    }

    // 'N' prior is the average:
    prior_set& nps(_lnprior[N_BASE]);
    for(unsigned i(0);i<N_BASE;++i) {
        prior_set& ps(_lnprior[i]);
        sum_gt(nps.genome,ps.genome);
        sum_gt(nps.poly,ps.poly);
    }
    norm_gt(nps.genome);
    norm_gt(nps.poly);

    // take logs:
    for(unsigned i(0);i<(N_BASE+1);++i) {
        prior_set& ps(_lnprior[i]);
        log_gt(ps.genome);
        log_gt(ps.poly);
    }  
}



static
void
increment_het_ratio_lhood(const extended_pos_info& epi,
                          const blt_float_t het_ratio,
                          blt_float_t* all_het_lhood,
                          const bool is_strand_specific,
                          const bool is_ss_fwd) {

    const blt_float_t chet_ratio(1.-het_ratio);

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

    const snp_pos_info& pi(epi.pi);
    const unsigned ref_gt(base_to_id(pi.ref_base));

    const unsigned n_calls(pi.calls.size());

    blt_float_t val_high[3];

    for(unsigned i(0);i<n_calls;++i){
        const blt_float_t eprob(epi.de[i]);
        const blt_float_t ceprob(1.-pi.calls[i].error_prob());

        // precalculate the result for expect values of 0.0, het_ratio, chet_ratio, 1.0
        val_high[0] = std::log(eprob)+log_one_third;
        val_high[1] = std::log((ceprob)*het_ratio+((1.-ceprob)*one_third)*chet_ratio);
        val_high[2] = std::log((ceprob)*chet_ratio+((1.-ceprob)*one_third)*het_ratio);

        const bool is_force_ref(is_strand_specific && (is_ss_fwd!=pi.calls[i].is_fwd_strand));

        const uint8_t obs_id(pi.calls[i].base_id);
        for(unsigned gt(N_BASE);gt<DIGT::SIZE;++gt){
            static const uint8_t low_remap[] = {0,2,1};
            const unsigned key(DIGT::expect2_bias(obs_id,(is_force_ref ? ref_gt : gt)));
            lhood_high[gt] += val_high[key];
            lhood_low[gt] += val_high[low_remap[key]];
        }
    }

    for(unsigned gt(0);gt<DIGT::SIZE;++gt){
        if(! DIGT::is_het(gt)) continue;
        all_het_lhood[gt] = log_sum(all_het_lhood[gt],lhood_high[gt]);
        all_het_lhood[gt] = log_sum(all_het_lhood[gt],lhood_low[gt]);
    }
}



void
pprob_digt_caller::
get_diploid_gt_lhood(const blt_options& opt,
                     const extended_pos_info& epi,
                     const bool is_het_bias,
                     const blt_float_t het_bias,
                     blt_float_t* const lhood,
                     const bool is_strand_specific,
                     const bool is_ss_fwd) {

    // get likelihood of each genotype
    for(unsigned gt(0);gt<DIGT::SIZE;++gt) lhood[gt] = 0.;

    const snp_pos_info& pi(epi.pi);
    const unsigned ref_gt(base_to_id(pi.ref_base));

    const unsigned n_calls(pi.calls.size());
    for(unsigned i(0);i<n_calls;++i){
        const base_call& bc(pi.calls[i]);
        const blt_float_t eprob(epi.de[i]);
        const blt_float_t ceprob(1.-bc.error_prob());
        const blt_float_t lnce(bc.ln_comp_error_prob());

        // precalculate the result for expect values of 0.0, 0.5 & 1.0
        blt_float_t val[3];
        val[0] = std::log(eprob)+log_one_third;
        val[1] = std::log((ceprob)+((1.-ceprob)*one_third))+log_one_half;
        val[2] = lnce;

        const bool is_force_ref(is_strand_specific && (is_ss_fwd!=bc.is_fwd_strand));

        const uint8_t obs_id(bc.base_id);
        for(unsigned gt(0);gt<DIGT::SIZE;++gt){
            lhood[gt] += val[DIGT::expect2(obs_id,(is_force_ref ? ref_gt : gt))];
        }
    }

    if(is_het_bias) {
        // loop is currently setup to assume a uniform het ratio subgenotype prior
        const unsigned n_bias_steps(1+static_cast<unsigned>(het_bias/opt.het_bias_max_ratio_inc));
        const blt_float_t ratio_increment(het_bias/static_cast<blt_float_t>(n_bias_steps));
        for(unsigned i(0);i<n_bias_steps;++i) {
            const blt_float_t het_ratio(0.5+(i+1)*ratio_increment);
            increment_het_ratio_lhood(epi,het_ratio,lhood,is_strand_specific,is_ss_fwd);
        }
        
        const unsigned n_het_subgt(1+2*n_bias_steps);
        const blt_float_t subgt_log_prior(std::log(1./static_cast<blt_float_t>(n_het_subgt)));
        
        for(unsigned gt(0);gt<DIGT::SIZE;++gt){
            if(! DIGT::is_het(gt)) continue;
            lhood[gt] += subgt_log_prior;
        }
    }


}



typedef diploid_genotype::result_set result_set;



void
debug_dump_digt_lhood(const blt_float_t* lhood,
                      std::ostream& os) {

    blt_float_t pprob[DIGT::SIZE];
    for(unsigned gt(0);gt<DIGT::SIZE;++gt){
        pprob[gt] = lhood[gt];
    }

    unsigned max_gt(0);
    normalize_ln_distro(pprob,pprob+DIGT::SIZE,max_gt);

    os << std::setprecision(3) << std::fixed;
    for(unsigned gt(0);gt<DIGT::SIZE;++gt){
        os << DIGT::label(gt) << ": " << -std::log(pprob[gt]) << " ";
    }
    os.unsetf(std::ios::fixed);
}



void
pprob_digt_caller::
calculate_result_set(const blt_float_t* lhood,
                     const blt_float_t* lnprior,
                     const unsigned ref_gt,
                     result_set& rs) {

    // mult by prior distro to get unnormalized pprob:
    //
    for(unsigned gt(0);gt<DIGT::SIZE;++gt){
        rs.pprob[gt] = lhood[gt] + lnprior[gt];
    }

    normalize_ln_distro(rs.pprob,rs.pprob+DIGT::SIZE,rs.max_gt);

    rs.snp_qphred=error_prob_to_qphred(rs.pprob[ref_gt]);
    rs.max_gt_qphred=error_prob_to_qphred(prob_comp(rs.pprob,rs.pprob+DIGT::SIZE,rs.max_gt));
}



///
/// The likelihood of an observed 'column' of observations for any
/// genotype: P(obs=ACT|genotype=AT) is a sum over all possible
/// columns:
///
/// sum { c in {A,C,G,T}**3 } { P(obs=ACT|col=c) * P(col=c|genotype=AT) }
///
/// Enumerating all {A,C,G,T}**col_size possible column values is
/// inefficent, so we only consider those columns for which
/// P(col|genotype) is nonzero.  For homozygous genotypes this always
/// leaves a single column; for heterozygous genotypes this is a sum
/// over all valid columns:
///
/// P(obs=ACT|genotype=AT)=
/// sum { c in {A,T}**3 } { P(obs=ACT|col=c) * P(col=c|genotype=AT }
///
/// that is...
///
/// P(obs|AAA)*P(AAA|genotype)+
/// P(obs|AAT)*P(AAT|genotype)+...etc
///
/// the first term is a product based on error probabilities:
/// P(obs=ACT|AAA) = (1-P(e1))*(P(e2)/3)*(P(e3)/3)
///
/// the second term is a product of genotype base frequencies:
/// P(AAA|genotype=AT) = .5**3
///
/// The likelihood calculated in the function below is an equivelent
/// factorization of the calculation described above.
///
void
pprob_digt_caller::
position_snp_call_pprob_digt(const blt_options& opt,
                             const extended_pos_info& epi,
                             diploid_genotype& dgt,
                             const bool is_always_test) const {

    const snp_pos_info& pi(epi.pi);

    if(pi.ref_base=='N') return;

    dgt.ref_gt=base_to_id(pi.ref_base);

    // check that a non-reference call meeting quality criteria even exists:
    if(! is_always_test) {
        if(is_spi_allref(pi,dgt.ref_gt)) return;
    }

    // get likelihood of each genotype
    blt_float_t lhood[DIGT::SIZE];
    get_diploid_gt_lhood(opt,epi,opt.is_bsnp_diploid_het_bias,opt.bsnp_diploid_het_bias,lhood);

    // get genomic site results:
    calculate_result_set(lhood,lnprior_genomic(dgt.ref_gt),dgt.ref_gt,dgt.genome);

    // get polymorphic site results:
    calculate_result_set(lhood,lnprior_polymorphic(dgt.ref_gt),dgt.ref_gt,dgt.poly);

    dgt.is_snp=(dgt.genome.snp_qphred != 0);

    // compute strand-bias here: 
    const bool is_compute_sb(true);
    if(is_compute_sb && dgt.is_snp) {
        blt_float_t lhood_fwd[DIGT::SIZE];
        get_diploid_gt_lhood(opt,epi,opt.is_bsnp_diploid_het_bias,opt.bsnp_diploid_het_bias,lhood_fwd,true,true);
        blt_float_t lhood_rev[DIGT::SIZE];
        get_diploid_gt_lhood(opt,epi,opt.is_bsnp_diploid_het_bias,opt.bsnp_diploid_het_bias,lhood_rev,true,false);

        // If max_gt is equal to reference, then go ahead and use it
        // for consistency, even though this makes the SB value
        // useless:
        const unsigned tgt(dgt.genome.max_gt);
        dgt.sb=std::max(lhood_fwd[tgt],lhood_rev[tgt])-lhood[tgt];
    }
}




std::ostream& operator<<(std::ostream& os,const diploid_genotype& dgt) {

    const result_set& ge(dgt.genome);
    const result_set& po(dgt.poly);

    os << " Q(snp): " << ge.snp_qphred
       << " max_gt: " << DIGT::label(ge.max_gt)
       << " Q(max_gt): "  << ge.max_gt_qphred
       << " max_gt|poly_site: " << po.max_gt
       << " Q(max_gt|poly_site): " << po.max_gt_qphred;

    return os;
}



void
write_diploid_genotype_allele(const blt_options& opt,
                              const snp_pos_info& pi,
                              const diploid_genotype& dgt,
                              std::ostream& os,
                              const bool is_nf_snp,
                              const double sb,
                              const unsigned hpol){

    const result_set& ge(dgt.genome);
    const result_set& po(dgt.poly);

    os << std::setprecision(6) << std::fixed;

    os << ge.snp_qphred
       << '\t' << DIGT::label(ge.max_gt) 
       << '\t' << ge.max_gt_qphred
       << '\t' << DIGT::label(po.max_gt)
       << '\t' << po.max_gt_qphred;

    if(opt.is_print_used_allele_counts) {
        pi.print_known_counts(os,opt.used_allele_count_min_qscore);
    }

    if(is_nf_snp) {
        os << '\t' << sb
           << '\t' << hpol;
    }

    if(opt.is_print_all_poly_gt) {
         for(unsigned gt(0);gt<DIGT::SIZE;++gt){
#if 1
             // print GT as prob:
             os << '\t' << po.pprob[gt];
#else
             // print GT as qval:
             os << '\t' << error_prob_to_qphred(prob_comp(po.pprob,po.pprob+DIGT::SIZE,gt));
#endif
         }
     }

     os.unsetf(std::ios::fixed);
}
