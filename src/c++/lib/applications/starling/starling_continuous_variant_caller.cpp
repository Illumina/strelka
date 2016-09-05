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



/// strand bias values summed over all samples
struct StrandBiasCounts
{
    unsigned fwdAlt = 0;
    unsigned revAlt = 0;
    unsigned fwdOther = 0;
    unsigned revOther = 0;
};



void
starling_continuous_variant_caller::
position_snp_call_continuous(
   const starling_base_options& opt,
   const snp_pos_info& good_pi,
   const unsigned baseId,
   GermlineContinuousSiteLocusInfo& locus)
{
    const uint8_t refBaseId = base_to_id(locus.ref);
    const bool isRefAllele(refBaseId == baseId);

    auto& siteAlleles(locus.getSiteAlleles());
    const auto alleleCount(siteAlleles.size());

    std::vector<StrandBiasCounts> strandBiasCounts(alleleCount);

    const unsigned sampleCount(locus.getSampleCount());
    for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
    {
        auto& sampleInfo(locus.getSample(sampleIndex));
        const auto& continuousSiteSampleInfo(locus.getContinuousSiteSample(sampleIndex));

        const double alleleFrequency(continuousSiteSampleInfo.getContinuousAlleleFrequency());

        {
            // use diploid gt codes as a convenient way to summarize the continuous variant calls:
            static const unsigned homrefGtIndex(VcfGenotypeUtil::getGenotypeIndex(0, 0));
            static const unsigned hetGtIndex(VcfGenotypeUtil::getGenotypeIndex(0, 1));
            static const unsigned homGtIndex(VcfGenotypeUtil::getGenotypeIndex(1, 1));

            auto getGtIndex = [&]() -> unsigned {
                if (isRefAllele)
                    return homrefGtIndex;
                else if (alleleFrequency >= (1. - opt.min_het_vf))
                    return homGtIndex;
                else if (alleleFrequency < opt.min_het_vf)
                    return homrefGtIndex; // STAR-66 - desired behavior
                else
                    return hetGtIndex;
            };

            sampleInfo.maxGenotypeIndexPolymorphic = getGtIndex();
        }

        sampleInfo.gqx = sampleInfo.genotypeQualityPolymorphic =
            poisson_qscore(continuousSiteSampleInfo.continuousAlleleDepth,
                           continuousSiteSampleInfo.continuousTotalDepth,
                           (unsigned) opt.min_qscore, 40);

        if (not isRefAllele)
        {
            // flag the whole site as a SNP if any call above the VF threshold is non-ref
            locus._is_snp = locus._is_snp || alleleFrequency > opt.min_het_vf;
        }

        // update strand bias intermediates:
        for (unsigned alleleIndex(0); alleleIndex < alleleCount; ++alleleIndex)
        {
            const auto& allele(siteAlleles[alleleIndex]);
            auto& sbcounts(strandBiasCounts[alleleIndex]);

            for (const base_call& bc : good_pi.calls)
            {
                if (bc.is_fwd_strand)
                {
                    if (bc.base_id == allele.baseId)
                        sbcounts.fwdAlt++;
                    else
                        sbcounts.fwdOther++;
                }
                else if (bc.base_id == allele.baseId)
                    sbcounts.revAlt++;
                else
                    sbcounts.revOther++;
            }
        }
    }

    // update strand bias for each allele:
    for (unsigned alleleIndex(0); alleleIndex < alleleCount; ++alleleIndex)
    {
        auto& allele(siteAlleles[alleleIndex]);
        const auto& sbcounts(strandBiasCounts[alleleIndex]);
        allele.strandBias = strand_bias(
            sbcounts.fwdAlt, sbcounts.revAlt, sbcounts.fwdOther, sbcounts.revOther, opt.noise_floor);
    }

    // get the qual score:
    locus.anyVariantAlleleQuality = 0;
    for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
    {
        const auto& sampleInfo(locus.getSample(sampleIndex));
        locus.anyVariantAlleleQuality = std::max(locus.anyVariantAlleleQuality, sampleInfo.genotypeQualityPolymorphic);
    }
}
