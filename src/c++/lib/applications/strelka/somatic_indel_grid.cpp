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

/// \file
/// \author Chris Saunders
///

#include "somatic_call_shared.hh"
#include "somatic_indel_grid.hh"
#include "qscore_calculator.hh"

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
somatic_indel_caller_grid::
somatic_indel_caller_grid(const strelka_options& opt)
    : _ln_som_match(log1p_switch(-opt.somatic_indel_rate)),
      _ln_som_mismatch(std::log(opt.somatic_indel_rate))
{
    calculateGermlineGenotypeLogPrior(opt.bindel_diploid_theta, _germlineGenotypeLogPrior);
}

static
void
get_indel_het_grid_lhood(const starling_base_options& opt,
                         const starling_base_deriv_options& dopt,
                         const starling_sample_options& sample_opt,
                         const IndelKey& indelKey,
                         const IndelSampleData& indelSampleData,
                         const bool is_include_tier2,
                         const bool is_use_alt_indel,
                         double* const lhood)
{
    static const unsigned lsize(DIGT_GRID::HET_RES*2);
    for (unsigned gt(0); gt<(lsize); ++gt) lhood[gt] = 0.;

    for (unsigned i(0); i<DIGT_GRID::HET_RES; ++i)
    {
        const double het_ratio((i+1)*DIGT_GRID::RATIO_INCREMENT);
        get_high_low_het_ratio_lhood(opt,dopt,
                                     sample_opt,
                                     indelKey,indelSampleData,het_ratio,
                                     is_include_tier2,is_use_alt_indel,
                                     lhood[lsize-(i+1)],
                                     lhood[i]);
    }
}



/// Test if the current target indel should be filtered because of other indels overlapping it.
///
/// This function will return true if the indel is not one of the top two indels by read support at the locus,
/// or if the top two alleles do not have 90% support among the top 3 overlapping alleles.
///
/// \return True if the indel should be filtered
static
bool
is_multi_indel_allele(
    const starling_base_deriv_options& dopt,
    const IndelSampleData& normalIndelSampleData,
    const IndelSampleData& tumorIndelSampleData,
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
    ReadPathScores total_pprob;
    get_sum_path_pprob(dopt,normalIndelSampleData,is_include_tier2,is_use_alt_indel,total_pprob,true);
    get_sum_path_pprob(dopt,tumorIndelSampleData,is_include_tier2,is_use_alt_indel,total_pprob,false);

    // next determine the top two indel alleles:
    std::vector<std::pair<double,int> > scores;
    scores.push_back(std::make_pair(-total_pprob.indel,static_cast<int>(INDEL)));
    scores.push_back(std::make_pair(-total_pprob.ref,static_cast<int>(REF)));
    const ReadPathScores::alt_indel_t& ai(total_pprob.alt_indel);
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
    // quickly find/filter the non-rejected overlapping indels:
    //
    is_overlap=((scores[0].second!=REF) && (scores[1].second!=REF));

    return false;
}



void
somatic_indel_caller_grid::
get_somatic_indel(
    const strelka_options& opt,
    const strelka_deriv_options& dopt,
    const starling_sample_options& normal_opt,
    const starling_sample_options& tumor_opt,
    const IndelKey& indelKey,
    const IndelData& indelData,
    const unsigned normalSampleIndex,
    const unsigned tumorSampleIndex,
    const bool is_use_alt_indel,
    somatic_indel_call& sindel) const
{
    // for now, lhood calculation of tumor and normal are independent:

    // get likelihood of each genotype
    double normal_lhood[DIGT_GRID::PRESTRAND_SIZE];
    double tumor_lhood[DIGT_GRID::PRESTRAND_SIZE];

    sindel.is_forced_output=indelData.isForcedOutput;

    const IndelSampleData& normalIndelSampleData(indelData.getSampleData(normalSampleIndex));
    const IndelSampleData& tumorIndelSampleData(indelData.getSampleData(tumorSampleIndex));

    static const unsigned n_tier(2);
    std::array<indel_result_set,n_tier> tier_rs;
    for (unsigned i(0); i<n_tier; ++i)
    {
        const bool is_include_tier2(i==1);
        if (is_include_tier2)
        {
            if (! opt.useTier2Evidence) continue;
            if (tier_rs[0].qphred==0)
            {
                if (! sindel.is_forced_output)   // if forced output then there's still a point to computing tier2
                {
                    tier_rs[1].qphred=0;
                    continue;
                }
            }
        }

#if 0
        std::cerr << "BUG: testing tier/ik: " << i << " " << ik;
#endif
        static const bool is_somatic_multi_indel_filter(true);
        bool isMultiIndelFilter(false);
        if (is_somatic_multi_indel_filter)
        {
            isMultiIndelFilter = (is_multi_indel_allele(dopt,normalIndelSampleData,tumorIndelSampleData,is_include_tier2,tier_rs[i].is_overlap));
            if (isMultiIndelFilter)
            {
#if 0
                std::cerr << "BUG: rejected\n";
#endif
                if (not sindel.is_forced_output)   // if forced output then there's still a point to computing tier2
                {
                    tier_rs[i].qphred = 0;
                    continue;
                }
            }
        }

        get_indel_digt_lhood(opt,dopt,normal_opt,
                             indelKey,normalIndelSampleData,
                             is_include_tier2,is_use_alt_indel,
                             normal_lhood);
        get_indel_digt_lhood(opt,dopt,tumor_opt,
                             indelKey,tumorIndelSampleData,
                             is_include_tier2,is_use_alt_indel,
                             tumor_lhood);

        get_indel_het_grid_lhood(opt,dopt,normal_opt,
                                 indelKey,normalIndelSampleData,
                                 is_include_tier2,is_use_alt_indel,
                                 normal_lhood+SOMATIC_DIGT::SIZE);
        get_indel_het_grid_lhood(opt,dopt,tumor_opt,
                                 indelKey,tumorIndelSampleData,
                                 is_include_tier2,is_use_alt_indel,
                                 tumor_lhood+SOMATIC_DIGT::SIZE);

        // TODO: this is a temporary solution
        blt_float_t normal_lhood_float[DIGT_GRID::PRESTRAND_SIZE];
        blt_float_t tumor_lhood_float[DIGT_GRID::PRESTRAND_SIZE];

        for (int j(0); j<DIGT_GRID::PRESTRAND_SIZE; ++j)
        {
            normal_lhood_float[j] = (blt_float_t) normal_lhood[j];
            tumor_lhood_float[j] = (blt_float_t) tumor_lhood[j];
        }

        // The somatic caller has not been configured to use sample-specific error rates yet, so arbitrarily pull
        // rates from the tumor-sample only for now.
        const double sharedIndelErrorRate(std::pow(tumorIndelSampleData.getErrorRates().indelToRefErrorProb.getValue(), opt.shared_indel_error_factor));
        const double logSharedIndelErrorRate(std::log(sharedIndelErrorRate)); // shared indel error rate
        const double logSharedIndelErrorRateComplement(log1p_switch(-sharedIndelErrorRate));

        calculate_result_set_grid(
            (float)opt.indel_contam_tolerance,
            (float)logSharedIndelErrorRate,
            (float)logSharedIndelErrorRateComplement,
            normal_lhood_float,
            tumor_lhood_float,
            _germlineGenotypeLogPrior,
            _ln_som_match,
            _ln_som_mismatch,
            tier_rs[i]
        );

        if (isMultiIndelFilter)
        {
            // this takes care of filtered cases which are also forced calls
            tier_rs[i].qphred = 0;
        }
    }

    if (! sindel.is_forced_output)
    {
        if (tier_rs[0].qphred==0 ||
            tier_rs[1].qphred==0) return;
    }

    sindel.sindel_tier=0;
    if (opt.useTier2Evidence)
    {
        if (tier_rs[0].qphred > tier_rs[1].qphred)
        {
            sindel.sindel_tier=1;
        }
    }

    sindel.sindel_from_ntype_tier=0;
    if (opt.useTier2Evidence)
    {
        if (tier_rs[0].from_ntype_qphred > tier_rs[1].from_ntype_qphred)
        {
            sindel.sindel_from_ntype_tier=1;
        }
    }

    sindel.rs=tier_rs[sindel.sindel_from_ntype_tier];

    if (tier_rs[0].ntype != tier_rs[1].ntype)
    {
        // catch NTYPE conflict states:
        sindel.rs.ntype = NTYPE::CONFLICT;
        sindel.rs.from_ntype_qphred = 0;

    }
    else
    {
        // classify NTYPE:
        //

        // convert diploid genotype into more limited ntype set:
        //
        if       (sindel.rs.ntype==SOMATIC_DIGT::REF)
        {
            sindel.rs.ntype=NTYPE::REF;
        }
        else if (sindel.rs.ntype==SOMATIC_DIGT::HOM)
        {
            sindel.rs.ntype=NTYPE::HOM;
        }
        else
        {
            sindel.rs.ntype=NTYPE::HET;
#if 0
        }
        else if (sindel.rs.ntype==SOMATIC_DIGT::HET)
        {
            sindel.rs.ntype=NTYPE::HET;
        }
        else
        {
            assert(0);
#endif
        }
    }

    sindel.rs.qphred = tier_rs[sindel.sindel_tier].qphred;
}
