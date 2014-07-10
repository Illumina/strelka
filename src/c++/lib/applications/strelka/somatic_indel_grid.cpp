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
    static const double error_mod( -std::log(static_cast<double>(STAR_DIINDEL_GRID::HET_RES*2)) );

    _ln_som_match=(log1p_switch(-opt.somatic_indel_rate));
    _ln_som_mismatch=(std::log(opt.somatic_indel_rate/(static_cast<double>(STAR_DIINDEL_GRID::SIZE-1))));

    std::fill(_lnprior.normal.begin(),_lnprior.normal.end(),0);
    std::fill(_lnprior.normal_poly.begin(),_lnprior.normal_poly.end(),0);

    const double ln_sie_rate( std::log(opt.shared_indel_error_rate) );
    const double ln_csie_rate( log1p_switch(-opt.shared_indel_error_rate) );

    const double* normal_lnprior_genomic(in_caller.lnprior_genomic());
    const double* normal_lnprior_polymorphic(in_caller.lnprior_polymorphic());
    for (unsigned ngt(0); ngt<STAR_DIINDEL::SIZE; ++ngt)
    {
        _lnprior.normal[ngt] = (normal_lnprior_genomic[ngt]+ln_csie_rate);
        _lnprior.normal_poly[ngt] = (normal_lnprior_polymorphic[ngt]+ln_csie_rate);
        //         _lnprior_normal_nonoise[ngt] = normal_lnprior[ngt];
    }

    for (unsigned ngt(STAR_DIINDEL::SIZE); ngt<STAR_DIINDEL_GRID::SIZE; ++ngt)
    {
        _lnprior.normal[ngt] = ln_sie_rate+error_mod;
        _lnprior.normal_poly[ngt] = ln_sie_rate+error_mod;
    }

#ifdef SOMATIC_DEBUG
    check_ln_distro(_lnprior_normal,
                    _lnprior_normal+STAR_DIINDEL_GRID::SIZE,
                    "somatic indel normal prior");
#endif
}



typedef somatic_indel_call::result_set result_set;



static
void
get_indel_het_grid_lhood(const starling_options& opt,
                         const starling_deriv_options& dopt,
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
void
calculate_result_set(const strelka_options& /*opt*/,
                     const std::vector<blt_float_t>& normal_lnprior,
                     const std::vector<blt_float_t>& /*normal_lnprior_poly*/,
                     const double lnmatch,
                     const double lnmismatch,
                     const double* normal_lhood,
                     const double* tumor_lhood,
                     result_set& rs)
{

#ifdef SOMATIC_DEBUG
    std::vector<double> check_prior(DDIINDEL_GRID::SIZE);

    for (unsigned ngt(0); ngt<STAR_DIINDEL_GRID::SIZE; ++ngt)
    {
        const double base_prior(normal_lnprior[ngt]);
        for (unsigned tgt(0); tgt<STAR_DIINDEL_GRID::SIZE; ++tgt)
        {
            const unsigned dgt(DDIINDEL_GRID::get_state(ngt,tgt));
            check_prior[dgt] =
                base_prior+
                ((tgt==ngt) ? lnmatch : lnmismatch);
        }
    }

    check_ln_distro(check_prior.begin(),
                    check_prior.end(),
                    "somatic indel full prior");
#endif

    // get unnormalized posterior:
    std::vector<double> pprob(DDIINDEL_GRID::SIZE);

    for (unsigned ngt(0); ngt<STAR_DIINDEL_GRID::SIZE; ++ngt)
    {
        const double base_prior(normal_lnprior[ngt]);
        for (unsigned tgt(0); tgt<STAR_DIINDEL_GRID::SIZE; ++tgt)
        {
            const unsigned dgt(DDIINDEL_GRID::get_state(ngt,tgt));
            pprob[dgt] =
                normal_lhood[ngt]+
                tumor_lhood[tgt]+
                base_prior+
                ((tgt==ngt) ? lnmatch : lnmismatch);
        }
    }

    //opt_normalize_ln_distro(pprob.begin(),pprob.end(),DDIINDEL_GRID::is_nonsom.val.begin(),rs.max_gt);
    normalize_ln_distro(pprob.begin(),pprob.end(),rs.max_gt);

#ifdef DEBUG_INDEL_CALL
    log_os << "INDEL_CALL pprob(noindel),pprob(hom),pprob(het): " << pprob[STAR_DIINDEL::NOINDEL] << " " << pprob[STAR_DIINDEL::HOM] << " " << pprob[STAR_DIINDEL::HET] << "\n";
#endif
    double nonsomatic_sum(0);
    for (unsigned gt(0); gt<STAR_DIINDEL_GRID::SIZE; ++gt)
    {
        nonsomatic_sum += pprob[DDIINDEL_GRID::get_state(gt,gt)];
    }
    rs.sindel_qphred=error_prob_to_qphred(nonsomatic_sum);

    if (0 == rs.sindel_qphred) return;

    // reset max_gt to the most likely state excluding normal noise states:
    //
    rs.max_gt=0;
    for (unsigned sgt(0); sgt<STAR_DIINDEL::SIZE; ++sgt)
    {
        for (unsigned tgt(0); tgt<STAR_DIINDEL_GRID::SIZE; ++tgt)
        {
            const unsigned dstate(DDIINDEL_GRID::get_state(sgt,tgt));
            if (pprob[dstate] > pprob[rs.max_gt]) rs.max_gt=dstate;
        }
    }

#ifndef USE_POLYNORMAL_SINDEL
    // now compute the probability that the event is notsomatic or notfrom each of
    // the three reference states:
    //
    double min_not_somfrom_sum(0);
    for (unsigned sgt(0); sgt<STAR_DIINDEL::SIZE; ++sgt)
    {
        double not_somfrom_sum(nonsomatic_sum);
        for (unsigned ngt(0); ngt<STAR_DIINDEL_GRID::SIZE; ++ngt)
        {
            // skip this case because we want events where the normal
            // state is NOT sgt:
            if (sgt==ngt) continue;
            for (unsigned tgt(0); tgt<STAR_DIINDEL_GRID::SIZE; ++tgt)
            {
                // skip this case because we've already summed
                // nonsomatic prob:
                if (tgt==ngt) continue;
                not_somfrom_sum += pprob[DDIINDEL_GRID::get_state(ngt,tgt)];
            }
        }
        if ((sgt==0) || (not_somfrom_sum<min_not_somfrom_sum))
        {
            min_not_somfrom_sum=not_somfrom_sum;
            rs.sindel_from_ntype_qphred=error_prob_to_qphred(not_somfrom_sum);
            rs.ntype=sgt;
        }
    }

    // double not_somfromanyhom_sum(nonsomatic_sum);
    // for(unsigned ngt(0);ngt<STAR_DIINDEL::SIZE;++ngt){
    //     if(STAR_DIINDEL::HET != ngt) continue;
    //     for(unsigned tgt(0);tgt<STAR_DIINDEL::SIZE;++tgt){
    //         if(tgt==ngt) continue;
    //         not_somfromanyhom_sum += pprob[DDIINDEL::get_state(ngt,tgt)];
    //     }
    // }
    // rs.sindel_from_anyhom_qphred=error_prob_to_qphred(not_somfromanyhom_sum);

    //    rs.max_gt_qphred=error_prob_to_qphred(prob_comp(pprob.begin(),pprob.end(),rs.max_gt));

#else

    // Calculate normal distribution alone so that we can classify this call:
    //
    // Polymorphic prior is used because in this situation we want to
    // be conservative about the reference classification --
    // ie. conditioned on only looking at putative somatic sites, we
    // require evidence to show that the normal is in fact reference
    // and not simply an unsampled copy of the somatic variation.
    //
    std::vector<double> normal_pprob(STAR_DIINDEL_GRID::SIZE);
    for (unsigned ngt(0); ngt<STAR_DIINDEL_GRID::SIZE; ++ngt)
    {
        normal_pprob[ngt] = normal_lhood[ngt]+normal_lnprior_poly[ngt];
    }

    unsigned max_norm_gt(0);
    normalize_ln_distro(normal_pprob.begin(),normal_pprob.end(),max_norm_gt);

    // find the probability of comp(max_norm_gt):
    const double ngt_prob(prob_comp(normal_pprob.begin(),normal_pprob.end(),max_norm_gt));

    // (1-(1-a)(1-b)) -> a+b-(ab)
    double not_somfrom_sum(nonsomatic_sum+ngt_prob-(nonsomatic_sum*ngt_prob));

    rs.sindel_from_ntype_qphred=error_prob_to_qphred(not_somfrom_sum);
    rs.ntype=max_norm_gt;
#endif
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
is_multi_indel_allele(const starling_deriv_options& dopt,
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
            if (! opt.is_tier2()) continue;
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

        calculate_result_set(opt,
                             lnprior_genomic(),
                             lnprior_polymorphic(),
                             _ln_som_match,_ln_som_mismatch,
                             normal_lhood,tumor_lhood,tier_rs[i]);
    }

    if (! sindel.is_forced_output)
    {
        if (tier_rs[0].sindel_qphred==0 ||
            tier_rs[1].sindel_qphred==0) return;
    }

    sindel.sindel_tier=0;
    if (opt.is_tier2())
    {
        if (tier_rs[0].sindel_qphred > tier_rs[1].sindel_qphred)
        {
            sindel.sindel_tier=1;
        }
    }

    sindel.sindel_from_ntype_tier=0;
    if (opt.is_tier2())
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
