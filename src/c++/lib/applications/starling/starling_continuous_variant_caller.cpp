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

#include "starling_continuous_variant_caller.hh"
#include "blt_util/qscore.hh"

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/distributions/binomial.hpp>



static double Likelihood(unsigned coverage, unsigned observedCallCount, double expectedFrequency)
{
    if (observedCallCount == 0)
        return 0;

    return boost::math::pdf(boost::math::binomial(coverage, expectedFrequency), observedCallCount);
}



static double AssignPValue(unsigned observedCallCount, unsigned coverage, unsigned estimatedBaseCallQuality)
{
    if (observedCallCount == 0)
        return 1.0;

    const double errorRate = qphred_to_error_prob(estimatedBaseCallQuality);
    return (boost::math::gamma_p(observedCallCount, coverage * errorRate));
}



int
starling_continuous_variant_caller::
poisson_qscore(unsigned callCount, unsigned coverage, unsigned estimatedBaseCallQuality, int maxQScore)
{
    double pValue = AssignPValue(callCount, coverage, estimatedBaseCallQuality);
    if (pValue <= 0) return maxQScore;
    return std::min(maxQScore, error_prob_to_qphred(pValue));
}




// calculate the ratio of the log likelihood of the variants on either strand / both strands
double starling_continuous_variant_caller::strand_bias(
    unsigned fwdAlt, unsigned revAlt, unsigned fwdOther, unsigned revOther, double /*noise*/)
{
    double expectedVf = (fwdAlt + revAlt) / ((double)fwdOther+revOther+fwdAlt+revAlt);

    // TODO: removed noise terms as they always evaluate to -inf
    auto fwd = std::log(Likelihood(fwdAlt+fwdOther, fwdAlt, expectedVf));// + std::log(Likelihood(revAlt+revOther, revAlt, noise));
    auto rev = std::log(Likelihood(revAlt+revOther, revAlt, expectedVf));// + std::log(Likelihood(fwdAlt+fwdOther, fwdAlt, noise));
    auto both = std::log(Likelihood(fwdAlt+fwdOther+revAlt+revOther, fwdAlt+revAlt, expectedVf));
    return std::max(fwd, rev) - both;
}

void starling_continuous_variant_caller::position_snp_call_continuous(
    const starling_base_options& opt,
    const snp_pos_info& good_pi,
    GermlineContinuousSiteLocusInfo& info)
{
    unsigned totalDepth = info.spanning_deletions;
    for (unsigned base_id(0); base_id<N_BASE; ++base_id)
    {
        totalDepth += info.alleleObservationCounts(base_id);
    }
    uint8_t ref_base_id = base_to_id(info.ref);

    auto generateCallInfo = [&](uint8_t base_id, bool force)
    {
        auto vf = info.alleleObservationCounts(base_id) / (double)totalDepth;
        if (vf > opt.min_het_vf || force)
        {
            GermlineContinuousSiteAlleleInfo call(totalDepth, info.alleleObservationCounts(base_id), (BASE_ID::index_t)base_id);
            call.gqx = call.gq = poisson_qscore(info.alleleObservationCounts(base_id), totalDepth, (unsigned)opt.min_qscore, 40);

            if (ref_base_id != base_id)
            {
                // flag the whole site as a SNP if any call above the VF threshold is non-ref
                info._is_snp = info._is_snp  || vf > opt.min_het_vf;
                unsigned int fwdAlt = 0;
                unsigned  revAlt = 0;
                unsigned  fwdOther = 0;
                unsigned  revOther = 0;
                for (const base_call& bc : good_pi.calls)
                {
                    if (bc.is_fwd_strand)
                    {
                        if (bc.base_id == base_id)
                            fwdAlt++;
                        else
                            fwdOther++;
                    }
                    else if (bc.base_id == base_id)
                        revAlt++;
                    else
                        revOther++;
                }

                call.strand_bias = strand_bias(fwdAlt, revAlt, fwdOther, revOther, opt.noise_floor);
            }
            info.altAlleles.push_back(call);
        }
    };


    for (unsigned base_id(0); base_id<N_BASE; ++base_id)
    {
        generateCallInfo(base_id, info.isForcedOutput);
    }
    if (info.altAlleles.empty())
    {
        // force at least a call for the reference so that we can assign filters to the locus (filters are in the calls)
        generateCallInfo(ref_base_id, true);
    }
}



void starling_continuous_variant_caller::add_indel_call(
    const starling_base_options& opt,
    const IndelKey& indelKey,
    const IndelData& indelData,
    const starling_indel_report_info& indelReportInfo,
    const starling_indel_sample_report_info& indelSampleReportInfo,
    GermlineContinuousIndelLocusInfo& info)
{
    // determine VF
    double vf = indelSampleReportInfo.n_confident_indel_reads / ((double)indelSampleReportInfo.total_confident_reads());
    if (vf > opt.min_het_vf || indelData.isForcedOutput)
    {
        info.altAlleles.emplace_back(
            indelSampleReportInfo.total_confident_reads(), indelSampleReportInfo.n_confident_indel_reads,
            indelKey, indelData, indelReportInfo, indelSampleReportInfo);
        GermlineContinuousIndelAlleleInfo& call = info.altAlleles.back();
        call.gqx = call.gq = poisson_qscore(indelSampleReportInfo.n_confident_indel_reads, indelSampleReportInfo.total_confident_reads(), (unsigned)opt.min_qscore, 40);
    }
    if (!info.altAlleles.empty())
    {
        info.is_het = info.altAlleles.size() > 1 || info.altAlleles.front().variant_frequency() < (1-opt.min_het_vf);
    }
}

