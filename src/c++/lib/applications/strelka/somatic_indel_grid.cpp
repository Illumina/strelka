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
#include "blt_util/digt.hh"
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
// ln_som_mismatch = log( P(S) )
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
    _ln_som_match=(log1p_switch(-opt.somatic_indel_rate));
    _ln_som_mismatch=(std::log(opt.somatic_indel_rate));

    std::fill(_bare_lnprior.normal.begin(),_bare_lnprior.normal.end(),0);

    const double* normal_lnprior_genomic(in_caller.lnprior_genomic());

    for (unsigned ngt(0); ngt<STAR_DIINDEL::SIZE; ++ngt)
    {
        _bare_lnprior.normal[ngt] = normal_lnprior_genomic[ngt];
    }
}

// Currently the index is reversed (i.e. index == 3 -> fraction == 0.95)
float get_fraction_from_index(int index)
{
    const float ratio_increment(0.5/static_cast<float>(STAR_DIINDEL_GRID::HET_RES+1));
    if(index == 0) return 0.f;
    if(index == 1) return 1.f;
    if(index == 2) return 0.5f;
    if(index < STAR_DIINDEL::SIZE+STAR_DIINDEL_GRID::HET_RES) return 1.0f - ratio_increment*(index-STAR_DIINDEL::SIZE+1);
    return 1.0f - ratio_increment*(index-STAR_DIINDEL::SIZE+2);
}

typedef somatic_indel_call::result_set result_set;

static
void
get_indel_het_grid_lhood(const starling_base_options& opt,
                         const starling_base_deriv_options& dopt,
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

void somatic_indel_caller_grid::set_somatic_prior(
        std::vector<double>& somatic_prior,
        const double ref_err_prob,
        const strelka_options& opt) const
{
    bool is_normal_contaminated = true;

    const double sie_rate(std::pow(ref_err_prob, opt.shared_indel_error_factor));
    const double ln_sie_rate(std::log(sie_rate));
    const double ln_csie_rate(log1p_switch(-sie_rate));

    double somatic_prior_normal[STAR_DIINDEL_GRID::SIZE] = {};
//    somatic_prior_normal[STAR_DIINDEL::NOINDEL] = 0.865894; // fn = 0.0
//    somatic_prior_normal[STAR_DIINDEL::HET] = 0.0;  // fn=0.5
//    somatic_prior_normal[STAR_DIINDEL::HOM] = 0.0;  // fn=1.0
//    somatic_prior_normal[STAR_DIINDEL_GRID::SIZE - 1] = 0.122517;  // fn = 0.05
//    somatic_prior_normal[STAR_DIINDEL_GRID::SIZE - 2] = 0.009934;  // fn = 0.1
//    somatic_prior_normal[STAR_DIINDEL_GRID::SIZE - 3] = 0.001656;  // fn = 0.15
    somatic_prior_normal[STAR_DIINDEL::NOINDEL] = 0.5; // fn = 0.0
    somatic_prior_normal[STAR_DIINDEL_GRID::SIZE - 1] = 0.5;  // fn = 0.05
//    somatic_prior_normal[STAR_DIINDEL_GRID::SIZE - 2] = 0.01;  // fn = 0.1

//    double somatic_prior_tumor[STAR_DIINDEL_GRID::SIZE] = {
//            0.001656, 0.0, 0.163907,    // 0.0, 1.0, 0.5
//            0.0, 0.0, 0.0, 0.0, 0.001656, 0.011589, 0.018212, 0.036424, 0.084437,    // 0.95, 0.90, ...
//            0.193709, 0.193709, 0.147351, 0.074503, 0.038079, 0.016556, 0.013245, 0.001656, 0.003311     // 0.45, 0.40, ...
//    };
    double somatic_prior_tumor[STAR_DIINDEL_GRID::SIZE];
    for(unsigned ft(0); ft<STAR_DIINDEL_GRID::SIZE; ++ft)
        somatic_prior_tumor[ft] = 1.0/static_cast<double>(STAR_DIINDEL_GRID::SIZE-1);

    for (unsigned ngt(0); ngt<STAR_DIINDEL::SIZE; ++ngt)
    {
        for (unsigned tgt(0); tgt<TWO_STATE_SOMATIC::SIZE; ++tgt) // 0: non-somatic, 1: somatic
        {
            for (unsigned ft(0); ft<STAR_DIINDEL_GRID::SIZE; ++ft)
            {
                for (unsigned fn(0); fn<STAR_DIINDEL_GRID::SIZE; ++fn)
                {
                    double lprob_f_given_g = 0.0;

                    if(tgt == 0)
                    {
                        if(fn == ft)
                        {
                            lprob_f_given_g = (fn == ngt) ? ln_csie_rate : ln_sie_rate-std::log(static_cast<double>(STAR_DIINDEL_GRID::SIZE-1));
                        }
                        else
                        {
                            lprob_f_given_g = -INFINITY;
                        }
                    }
                    else
                    {
                        if(fn == ft)
                        {
                            lprob_f_given_g = -INFINITY;
                        }
                        else
                        {
                            if(!is_normal_contaminated || ngt != STAR_DIINDEL::NOINDEL)
                            {
                                lprob_f_given_g = -std::log(static_cast<double>(STAR_DIINDEL_GRID::SIZE-1));
                                lprob_f_given_g += (fn == ngt) ? ln_csie_rate : ln_sie_rate-std::log(static_cast<double>(STAR_DIINDEL_GRID::SIZE-1));
                            }
                            else
                            {
                                // fn should be smaller than ft
                                if(get_fraction_from_index(fn) >= get_fraction_from_index(ft))
                                    lprob_f_given_g = -INFINITY;
                                else
                                {
                                    lprob_f_given_g = std::log(somatic_prior_normal[fn]) + std::log(somatic_prior_tumor[ft]);
                                }
                            }
                        }
                    }

                    const unsigned dgt(DDIINDEL_GRID::get_state(fn, ft));
                    const unsigned genotype_index(tgt*STAR_DIINDEL::SIZE + ngt);
                    const unsigned index(genotype_index*DDIINDEL_GRID::SIZE + dgt);
                    somatic_prior[index] = lprob_f_given_g;
                }
            }
        }
    }
}

void
somatic_indel_caller_grid::calculate_result_set(
    const std::vector<double>& somatic_prior,
    const double* normal_lhood,
    const double* tumor_lhood,
    result_set& rs) const
{
    double log_post_prob[STAR_DIINDEL::SIZE][TWO_STATE_SOMATIC::SIZE];
    double max_log_prob = -INFINITY;
    rs.max_gt=0;

    for (unsigned ngt(0); ngt<STAR_DIINDEL::SIZE; ++ngt)
    {
        // logP(Gn=ngt, Gt=tgt)
        double log_diploid_prior_prob = _bare_lnprior.normal[ngt];  // logP(Gn)
        for (unsigned tgt(0); tgt<TWO_STATE_SOMATIC::SIZE; ++tgt) // 0: non-somatic, 1: somatic
        {
            double log_prior_prob = log_diploid_prior_prob + ((tgt == 0) ? _ln_som_match : _ln_som_mismatch);
            double max_log_sum = -INFINITY;
            double log_sum[DDIINDEL_GRID::SIZE];
            for (unsigned ft(0); ft<STAR_DIINDEL_GRID::SIZE; ++ft)
            {
                for (unsigned fn(0); fn<STAR_DIINDEL_GRID::SIZE; ++fn)
                {
                    const unsigned dgt(DDIINDEL_GRID::get_state(fn, ft));
                    const unsigned genotype_index(tgt*STAR_DIINDEL::SIZE + ngt);
                    const unsigned index(genotype_index*DDIINDEL_GRID::SIZE + dgt);

                    double sum = somatic_prior[index] + normal_lhood[fn] + tumor_lhood[ft];
                    log_sum[dgt] = sum;

                    if(sum > max_log_sum) max_log_sum = sum;
                }
            }

            // calculate log(exp(log_sum[0])+exp(log_sum[1])+...)
            double sum = 0.0;
            for (int i(0); i<STAR_DIINDEL_GRID::SIZE*STAR_DIINDEL_GRID::SIZE; ++i)
            {
                sum += std::exp(log_sum[i] - max_log_sum);
            }

            log_post_prob[ngt][tgt] = log_prior_prob + max_log_sum + std::log(sum);


            if(log_post_prob[ngt][tgt] > max_log_prob)
            {
                max_log_prob = log_post_prob[ngt][tgt];
                rs.max_gt = DDIINDEL_GRID::get_state(ngt,tgt);  // TODO: this can be simplified
            }
        }
    }

    //#ifdef DEBUG_INDEL_CALL
    //    log_os << "INDEL_CALL pprob(noindel),pprob(hom),pprob(het): " << pprob[STAR_DIINDEL::NOINDEL] << " " << pprob[STAR_DIINDEL::HOM] << " " << pprob[STAR_DIINDEL::HET] << "\n";
    //#endif

    double post_prob[STAR_DIINDEL::SIZE][TWO_STATE_SOMATIC::SIZE];
    double sum_prob = 0.0;
    for (unsigned ngt(0); ngt<STAR_DIINDEL::SIZE; ++ngt)
    {
        for (unsigned tgt(0); tgt<TWO_STATE_SOMATIC::SIZE; ++tgt)
        {
            double prob = std::exp(log_post_prob[ngt][tgt] - max_log_prob); // to prevent underflow
            sum_prob += prob;
        }
    }

    // normalize probabilities, determine QSI, NT and QSI_NT
    // Compute QSI_NT
    //
    // compute the probability that the indel is:
    //
    // not_somfrom_sum =
    //     ((not somatic) OR (not 'from one of the three reference states'))
    //
    // ...where we minimize over the values from the 3 ref states
    //
    // We're interested in this value b/c it is the complement of:
    // (somatic AND 'from one of the three reference states')
    //
    // which is the value we want to provide a quality score for --
    // the quality score is an error term, so it is  based on the
    // the complement value.
    //
    //
    // NOTE this is unnecessarily inefficient, even after controlling
    // for potential round-off error (which itself might be overkill), we're
    // still repeatedly summing most of the frequency space 3 times.
    //
    double log_sum_prob = std::log(sum_prob);
    double min_not_somfrom_sum(INFINITY);
    double nonsom_prob = 0.0;
    for (unsigned ngt(0); ngt<STAR_DIINDEL::SIZE; ++ngt)
    {
        double som_prob_given_ngt(0);
        for (unsigned tgt(0); tgt<TWO_STATE_SOMATIC::SIZE; ++tgt)
        {
            post_prob[ngt][tgt] = std::exp(log_post_prob[ngt][tgt] - max_log_prob - log_sum_prob);
            if(tgt == 0)
            {
                nonsom_prob += post_prob[ngt][tgt];
            }
            else
            {
                som_prob_given_ngt += post_prob[ngt][tgt];
            }
        }

        double err_som_and_ngt = 1.0 - som_prob_given_ngt;
        if (err_som_and_ngt<min_not_somfrom_sum)
        {
            min_not_somfrom_sum=err_som_and_ngt;
            rs.sindel_from_ntype_qphred=error_prob_to_qphred(err_som_and_ngt);
            rs.ntype=ngt;
        }
    }
    rs.sindel_qphred=error_prob_to_qphred(nonsom_prob);

    if (0 == rs.sindel_qphred) return;
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
is_multi_indel_allele(const starling_base_deriv_options& dopt,
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
    std::vector<double> somatic_prior(TWO_STATE_SOMATIC::SIZE*STAR_DIINDEL::SIZE*DDIINDEL_GRID::SIZE, 0);

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
            if (! opt.tier2.is_tier2()) continue;
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
                                                normal_lhood
                                                );
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

        set_somatic_prior(somatic_prior, ref_error_prob, opt);

        calculate_result_set(somatic_prior,
                normal_lhood,tumor_lhood,tier_rs[i]);
    }

    if (! sindel.is_forced_output)
    {
        if (tier_rs[0].sindel_qphred==0 ||
            tier_rs[1].sindel_qphred==0) return;
    }

    sindel.sindel_tier=0;
    if (opt.tier2.is_tier2())
    {
        if (tier_rs[0].sindel_qphred > tier_rs[1].sindel_qphred)
        {
            sindel.sindel_tier=1;
        }
    }

    sindel.sindel_from_ntype_tier=0;
    if (opt.tier2.is_tier2())
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
