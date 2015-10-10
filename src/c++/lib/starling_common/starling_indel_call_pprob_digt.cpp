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

///
/// \author Chris Saunders
///

#include "blt_util/math_util.hh"
#include "blt_util/prob_util.hh"
#include "starling_common/indel_util.hh"
#include "starling_common/starling_indel_call_pprob_digt.hh"
#include "starling_common/starling_indel_report_info.hh"

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <iostream>
#include <iomanip>


//#define DEBUG_INDEL_CALL

#ifdef DEBUG_INDEL_CALL
#include "blt_util/log.hh"
#endif



static
void
set_diploid_prior(
    const double theta,
    indel_digt_caller::prior_group& prior)
{
    using namespace STAR_DIINDEL;

    prior.genome[NOINDEL]=log1p_switch(-(3.*theta)/2.);
    prior.genome[HOM]=std::log(theta/2.);
    prior.genome[HET]=std::log(theta);

    prior.poly[NOINDEL]=std::log(0.25);
    prior.poly[HOM]=std::log(0.25);
    prior.poly[HET]=std::log(0.5);
}



static
void
set_haploid_prior(
    const double theta,
    indel_digt_caller::prior_group& prior)
{
    static const double log0(-std::numeric_limits<double>::infinity());

    using namespace STAR_DIINDEL;

    prior.genome[NOINDEL]=log1p_switch(-theta);
    prior.genome[HOM]=std::log(theta);
    prior.genome[HET]=log0;

    prior.poly[NOINDEL]=std::log(0.5);
    prior.poly[HOM]=std::log(0.5);
    prior.poly[HET]=log0;
}



indel_digt_caller::
indel_digt_caller(const double theta)
{
    set_diploid_prior(theta,_lnprior);
    set_haploid_prior(theta,_lnprior_haploid);
}



static
double
integrate_out_sites(const starling_base_deriv_options& dopt,
                    const uint16_t nsite,
                    const double p_on_site,
                    const bool is_tier2_pass)
{
    return log_sum((p_on_site + dopt.site_lnprior),
                   (dopt.get_nonsite_path_lnp(is_tier2_pass,nsite) + dopt.nonsite_lnprior));
}



// Calculate the expected ratio of reads observed to support each
// allele. Note that for sites and single breakpoints this is expected
// to match the sample allele ratio, however for indels this can
// change as a funciton of indel and read length.
//
// Note this routine does not account for overlapping indels
//
static
void
get_het_observed_allele_ratio(const unsigned read_length,
                              const unsigned min_overlap,
                              const indel_key& ik,
                              const double het_allele_ratio,
                              double& log_ref_prob,
                              double& log_indel_prob)
{
    assert((ik.type==INDEL::INSERT) ||
           (ik.type==INDEL::DELETE) ||
           (ik.type == INDEL::SWAP));

    // the expected relative read depth for two breakpoints separated by a distance of 0:
    const unsigned base_expect( (read_length+1)<(2*min_overlap) ? 0 : (read_length+1)-(2*min_overlap) );

    // Get expected relative read depth for the shorter and longer
    // paths of a general sequence replacement. Note this includes
    // basic insertions and deletions, in these cases
    // spath_break_distance is 0 and spath_expect equals base_expect:
    //
    const double ref_path_expect(base_expect+std::min(ik.delete_length(),base_expect));
    const double indel_path_expect(base_expect+std::min(ik.insert_length(),base_expect));
    const double ref_path_term((1-het_allele_ratio)*ref_path_expect);
    const double indel_path_term(het_allele_ratio*indel_path_expect);
    const double total_path_term(ref_path_term+indel_path_term);

    if (total_path_term>0)
    {
        const double indel_prob(indel_path_term/total_path_term);
        log_ref_prob=std::log(1.-indel_prob);
        log_indel_prob=std::log(indel_prob);
    }
}



void
indel_digt_caller::
get_high_low_het_ratio_lhood(
    const starling_base_options& /*opt*/,
    const starling_base_deriv_options& dopt,
    const starling_sample_options& sample_opt,
    const double indel_error_lnp,
    const double indel_real_lnp,
    const double ref_error_lnp,
    const double ref_real_lnp,
    const indel_key& ik,
    const indel_data& id,
    const double het_ratio,
    const bool is_tier2_pass,
    const bool is_use_alt_indel,
    double& het_lhood_high,
    double& het_lhood_low)
{
    // handle het ratio and its complement in one step:
    const double chet_ratio(1.-het_ratio);

    const double log_het_ratio(std::log(het_ratio));
    const double log_chet_ratio(std::log(chet_ratio));

    const bool is_breakpoint(ik.is_breakpoint());

    het_lhood_high=0;
    het_lhood_low=0;

    for (const auto& score : id.read_path_lnp)
    {
        const read_path_scores& path_lnp(score.second);

        // optionally skip tier2 data:
        if ((! is_tier2_pass) && (! path_lnp.is_tier1_read)) continue;

        // get alt path lnp:
        double alt_path_lnp(path_lnp.ref);
#if 0
        if (is_use_alt_indel && path_lnp.is_alt &&
            (path_lnp.alt > alt_path_lnp))
        {
            alt_path_lnp=path_lnp.alt;
        }
#else
        if (is_use_alt_indel)
        {
            for (const auto& alt : path_lnp.alt_indel)
            {
                if (alt.second>alt_path_lnp) alt_path_lnp=alt.second;
            }
        }
#endif

        const double noindel_lnp(log_sum(alt_path_lnp+ref_real_lnp,path_lnp.indel+indel_error_lnp));
        const double hom_lnp(log_sum(alt_path_lnp+ref_error_lnp,path_lnp.indel+indel_real_lnp));

        // allele ratio convention is that the indel occurs at the
        // het_allele ratio and the alternate allele occurs at
        // (1-het_allele_ratio):
        {
            double log_ref_prob(log_chet_ratio);
            double log_indel_prob(log_het_ratio);
            if (! is_breakpoint)
            {
                get_het_observed_allele_ratio(path_lnp.read_length,sample_opt.min_read_bp_flank,
                                              ik,het_ratio,log_ref_prob,log_indel_prob);
            }
            const double het_lnp(log_sum(noindel_lnp+log_ref_prob,hom_lnp+log_indel_prob));

            het_lhood_low += integrate_out_sites(dopt,path_lnp.nsite,het_lnp,is_tier2_pass);
        }

        {
            double log_ref_prob(log_het_ratio);
            double log_indel_prob(log_chet_ratio);
            if (! is_breakpoint)
            {
                get_het_observed_allele_ratio(path_lnp.read_length,sample_opt.min_read_bp_flank,
                                              ik,chet_ratio,log_ref_prob,log_indel_prob);
            }
            const double het_lnp(log_sum(noindel_lnp+log_ref_prob,hom_lnp+log_indel_prob));

            het_lhood_high += integrate_out_sites(dopt,path_lnp.nsite,het_lnp,is_tier2_pass);
        }
    }
}



static
void
increment_het_ratio_lhood(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const starling_sample_options& sample_opt,
    const double indel_error_lnp,
    const double indel_real_lnp,
    const double ref_error_lnp,
    const double ref_real_lnp,
    const indel_key& ik,
    const indel_data& id,
    const double het_ratio,
    const bool is_tier2_pass,
    const bool is_use_alt_indel,
    double* const lhood)
{
    // high and low allele ratio variants:
    double het_lhood_high;
    double het_lhood_low;

    indel_digt_caller::get_high_low_het_ratio_lhood(opt,dopt,sample_opt,
                                                    indel_error_lnp,indel_real_lnp,
                                                    ref_error_lnp,ref_real_lnp,
                                                    ik,id,het_ratio,is_tier2_pass,
                                                    is_use_alt_indel,
                                                    het_lhood_high,het_lhood_low);

    lhood[STAR_DIINDEL::HET] = log_sum(lhood[STAR_DIINDEL::HET],het_lhood_low);
    lhood[STAR_DIINDEL::HET] = log_sum(lhood[STAR_DIINDEL::HET],het_lhood_high);
}



// total the path likelihoods of ref,indel and alt_indel states
//
void
get_sum_path_pprob(
    const starling_base_deriv_options& dopt,
    const indel_data& id,
    const bool is_tier2_pass,
    const bool is_use_alt_indel,
    read_path_scores& total_pprob,
    const bool is_init_total)
{
    static const double initval(0);

    if (is_init_total)
    {
        total_pprob.ref=initval;
        total_pprob.indel=initval;
        total_pprob.nsite=0;
    }

    typedef std::map<indel_key,unsigned> aimap_t;
    aimap_t alt_indel_index;

    for (const auto& score : id.read_path_lnp)
    {
        const read_path_scores& path_lnp(score.second);

        // optionally skip tier2 data:
        if ((! is_tier2_pass) && (! path_lnp.is_tier1_read)) continue;

        const read_path_scores path_pprob(indel_lnp_to_pprob(dopt,path_lnp,is_tier2_pass,is_use_alt_indel));

        total_pprob.indel += path_pprob.indel;
        total_pprob.ref += path_pprob.ref;

        if (! is_use_alt_indel) continue;

        for (const auto& alt : path_pprob.alt_indel)
        {
            aimap_t::iterator tj(alt_indel_index.find(alt.first));
            if (tj==alt_indel_index.end())
            {
                alt_indel_index[alt.first]=total_pprob.alt_indel.size();
                total_pprob.alt_indel.push_back(alt);
            }
            else
            {
                total_pprob.alt_indel[tj->second].second += alt.second;
            }
        }
    }
}



// Run a simple approx noise-filter first for sites with many indel
// alleles
//
// The goal of this filter is to identify obvious invalid cases where
// these are detectable without a full haplotype model -- it thus will
// let some obvious noise cases through.
//
// is_include_alt_indel == true is implied by calling this function
//
static
bool
is_diploid_indel_noise(
    const starling_base_deriv_options& dopt,
    const indel_data& id,
    const bool is_tier2_pass)
{
    static const bool is_use_alt_indel(true);

    // test is to *sum* the likelihoods supporting each indel:
    //

    // first sum every read intersecting the indel:
    read_path_scores total_pprob;
    get_sum_path_pprob(dopt,id,is_tier2_pass,is_use_alt_indel,total_pprob);

    enum iallele_t
    {
        INDEL = -2,
        REF = -1
    };

    // next determine the top two indel alleles:
    int max1_id(INDEL);
    int max2_id(REF);
    double max1(total_pprob.indel);
    double max2(total_pprob.ref);
    if (max1<max2)
    {
        std::swap(max1,max2);
        std::swap(max1_id,max2_id);
    }
    const read_path_scores::alt_indel_t& ai(total_pprob.alt_indel);
    const int ais(ai.size());
    for (int i(0); i<ais; ++i)
    {
        if       (ai[i].second>max1)
        {
            max2=max1;
            max2_id=max1_id;
            max1=ai[i].second;
            max1_id=i;
        }
        else if (ai[i].second>max2)
        {
            max2=ai[i].second;
            max2_id=i;
        }
    }

#if 0
    std::cerr << "BUG: id: " << id;
    std::cerr << "BUG: total_pprob: " << total_pprob << "\n";
    std::cerr << "BUG: max1_id,max2_id " << max1_id << " " << max2_id << "\n";
#endif

    // If the top two alleles are both alternate indels, check that
    // they interfere with each other.  If not, we are forced to make
    // the conservative assumption that they occur as part of the same
    // haplotype:
    //
    if (max1_id>=0 && max2_id>=0)
    {
        if (! is_indel_conflict(ai[max1_id].first,ai[max2_id].first))
        {
            return (total_pprob.ref>total_pprob.indel);
        }
    }

    return ((max1_id!=INDEL) && (max2_id!=INDEL));
}



void
indel_digt_caller::
get_indel_digt_lhood(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const starling_sample_options& sample_opt,
    const double indel_error_prob,
    const double ref_error_prob,
    const indel_key& ik,
    const indel_data& id,
    const bool is_het_bias,
    const double het_bias,
    const bool is_tier2_pass,
    const bool is_use_alt_indel,
    double* const lhood)
{
    static const double loghalf(-std::log(2.));

    for (unsigned gt(0); gt<STAR_DIINDEL::SIZE; ++gt) lhood[gt] = 0.;

    const bool is_breakpoint(ik.is_breakpoint());

    const double indel_error_lnp(std::log(indel_error_prob));
    const double indel_real_lnp(std::log(1.-indel_error_prob));
    const double ref_error_lnp(std::log(ref_error_prob));
    const double ref_real_lnp(std::log(1.-ref_error_prob));

    for (const auto& score : id.read_path_lnp)
    {
        const read_path_scores& path_lnp(score.second);

        // optionally skip tier2 data:
        if ((! is_tier2_pass) && (! path_lnp.is_tier1_read)) continue;

        // get alt path lnp:
        double alt_path_lnp(path_lnp.ref);
#if 0
        if (is_use_alt_indel && path_lnp.is_alt &&
            (path_lnp.alt > alt_path_lnp))
        {
            alt_path_lnp=path_lnp.alt;
        }
#else
        if (is_use_alt_indel)
        {
            for (const auto& alt : path_lnp.alt_indel)
            {
                if (alt.second>alt_path_lnp) alt_path_lnp=alt.second;
            }
        }
#endif

        const double noindel_lnp(log_sum(alt_path_lnp+ref_real_lnp,path_lnp.indel+indel_error_lnp));
        const double hom_lnp(log_sum(alt_path_lnp+ref_error_lnp,path_lnp.indel+indel_real_lnp));

        // allele ratio convention is that the indel occurs at the
        // het_allele ratio and the alternate allele occurs at
        // (1-het_allele_ratio):

        double log_ref_prob(loghalf);
        double log_indel_prob(loghalf);
        if (! is_breakpoint)
        {
            static const double het_allele_ratio(0.5);
            get_het_observed_allele_ratio(path_lnp.read_length,sample_opt.min_read_bp_flank,
                                          ik,het_allele_ratio,log_ref_prob,log_indel_prob);
        }
        const double het_lnp(log_sum(noindel_lnp+log_ref_prob,hom_lnp+log_indel_prob));

        lhood[STAR_DIINDEL::NOINDEL] += integrate_out_sites(dopt,path_lnp.nsite,noindel_lnp,is_tier2_pass);
        lhood[STAR_DIINDEL::HOM]     += integrate_out_sites(dopt,path_lnp.nsite,hom_lnp,is_tier2_pass);
        lhood[STAR_DIINDEL::HET]     += integrate_out_sites(dopt,path_lnp.nsite,het_lnp,is_tier2_pass);

#ifdef DEBUG_INDEL_CALL
        //log_os << std::setprecision(8);
        //log_os << "INDEL_CALL i,ref_lnp,indel_lnp,lhood(noindel),lhood(hom),lhood(het): " << i << " " << path_lnp.ref << " " << path_lnp.indel << " " << lhood[STAR_DIINDEL::NOINDEL] << " " << lhood[STAR_DIINDEL::HOM] << " " << lhood[STAR_DIINDEL::HET] << "\n";
#endif
    }


    if (is_het_bias)
    {
        // loop is currently setup to assume a uniform het ratio subgenotype prior
        const unsigned n_bias_steps(1+static_cast<unsigned>(het_bias/opt.het_bias_max_ratio_inc));
        const double ratio_increment(het_bias/static_cast<double>(n_bias_steps));
        for (unsigned step(0); step<n_bias_steps; ++step)
        {
            const double het_ratio(0.5+(step+1)*ratio_increment);
            increment_het_ratio_lhood(opt,dopt,sample_opt,
                                      indel_error_lnp,indel_real_lnp,
                                      ref_error_lnp,ref_real_lnp,
                                      ik,id,het_ratio,is_tier2_pass,is_use_alt_indel,lhood);
        }

        const unsigned n_het_subgt(1+2*n_bias_steps);
        const double subgt_log_prior(std::log(static_cast<double>(n_het_subgt)));
        lhood[STAR_DIINDEL::HET] -= subgt_log_prior;
    }
}



///
///
void
indel_digt_caller::
starling_indel_call_pprob_digt(const starling_base_options& opt,
                               const starling_base_deriv_options& dopt,
                               const starling_sample_options& sample_opt,
                               const double indel_error_prob,
                               const double ref_error_prob,
                               const indel_key& ik,
                               const indel_data& id,
                               const bool is_use_alt_indel,
                               starling_diploid_indel& dindel) const
{
    // no immediate plans to include this for regular indel-calling:
    static const bool is_tier2_pass(false);

    const bool is_haploid(dindel.is_haploid());

    if (opt.is_noise_indel_filter && is_diploid_indel_noise(dopt,id,is_tier2_pass))
    {
        dindel.is_indel=false;
        return;
    }

    // turn off het bias in haploid case:
    const bool is_het_bias((!is_haploid) && opt.is_bindel_diploid_het_bias);

    // get likelihood of each genotype:
    double lhood[STAR_DIINDEL::SIZE];
    get_indel_digt_lhood(opt,dopt,sample_opt,indel_error_prob,ref_error_prob,ik,id,
                         is_het_bias,opt.bindel_diploid_het_bias,
                         is_tier2_pass,is_use_alt_indel,lhood);

    // mult by prior distro to get unnormalized pprob:
    {
        const double* indel_lnprior(lnprior_genomic(is_haploid));
        for (unsigned gt(0); gt<STAR_DIINDEL::SIZE; ++gt)
        {
            dindel.pprob[gt] = lhood[gt] + indel_lnprior[gt];
        }
    }

    normalize_ln_distro(dindel.pprob,dindel.pprob+STAR_DIINDEL::SIZE,dindel.max_gt);

#ifdef DEBUG_INDEL_CALL
    log_os << "INDEL_CALL pprob(noindel),pprob(hom),pprob(het): " << dindel.pprob[STAR_DIINDEL::NOINDEL] << " " << dindel.pprob[STAR_DIINDEL::HOM] << " " << dindel.pprob[STAR_DIINDEL::HET] << "\n";
#endif

    dindel.indel_qphred=error_prob_to_qphred(dindel.pprob[STAR_DIINDEL::NOINDEL]);
    dindel.max_gt_qphred=error_prob_to_qphred(prob_comp(dindel.pprob,dindel.pprob+STAR_DIINDEL::SIZE,dindel.max_gt));

    // set phredLoghood:
    {
        unsigned gtcount(STAR_DIINDEL::SIZE);
        if (is_haploid) gtcount=2;
        unsigned maxIndex(0);
        for (unsigned gt(1); gt<gtcount; ++gt)
        {
            if (lhood[gt] > lhood[maxIndex]) maxIndex = gt;
        }
        for (unsigned gt(0); gt<gtcount; ++gt)
        {
            dindel.phredLoghood[gt] = std::min(dindel.maxQ,ln_error_prob_to_qphred(lhood[gt]-lhood[maxIndex]));
        }
    }

    // add new poly calls (this trashes lhood):
    {
        const double* indel_lnprior(lnprior_polymorphic(is_haploid));
        for (unsigned gt(0); gt<STAR_DIINDEL::SIZE; ++gt)
        {
            lhood[gt] += indel_lnprior[gt];
        }
    }
    normalize_ln_distro(lhood,lhood+STAR_DIINDEL::SIZE,dindel.max_gt_poly);
    dindel.max_gt_poly_qphred=error_prob_to_qphred(prob_comp(lhood,lhood+STAR_DIINDEL::SIZE,dindel.max_gt_poly));

    // old report criteria:
    //dindel.is_indel=(dindel.max_gt != STAR_DIINDEL::NOINDEL);
    dindel.is_indel=(dindel.indel_qphred != 0);
}
