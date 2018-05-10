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


#include "gvcf_locus_info.hh"
#include "blt_util/math_util.hh"
#include "common/Exceptions.hh"

#include "boost/math/distributions/binomial.hpp"
#include "boost/math/distributions.hpp"
#include "rnaVariantEmpiricalScoringFeatures.hh"
#include "starling_continuous_variant_caller.hh"

#include <iostream>
#include <map>
#include <sstream>
#include <typeinfo>



void
GermlineFilterKeeper::
write(std::ostream& os) const
{
    if (filters.none())
    {
        os << "PASS";
        return;
    }

    bool is_sep(false);
    for (unsigned i(0); i<GERMLINE_VARIANT_VCF_FILTERS::SIZE; ++i)
    {
        if (! filters.test(i)) continue;

        if (is_sep)
        {
            os << ";";
        }
        else
        {
            is_sep=true;
        }
        os << GERMLINE_VARIANT_VCF_FILTERS::get_label(i);
    }
}



std::ostream&
operator<<(
    std::ostream& os,
    const LocusSampleInfo& lsi)
{
    os << "LocusSampleInfo ploidy: " << lsi.getPloidy()
       << " maxGt: " << lsi.maxGenotypeIndexPolymorphic
       << " gq: " << lsi.genotypeQualityPolymorphic
       << " gqx: " << lsi.gqx
       << " filters: ";
    lsi.filters.write(os);
    os << "\n";
    os << "supportCounts: " << lsi.supportCounts;
    os << "PL: ";
    bool isFirst(true);
    for (const auto pl : lsi.genotypePhredLoghood)
    {
        if (not isFirst) os << ':';
        os << pl;
        isFirst = false;
    }
    os << '\n';
    return os;
}



std::ostream&
operator<<(
    std::ostream& os,
    const LocusInfo& li)
{
    os << "LocusInfo pos: " << li.pos
       << " qual: " << li.anyVariantAlleleQuality;

    const unsigned sampleCount(li.getSampleCount());
    os << " sampleCount: " << sampleCount << "\n";

    for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
    {
        os << "SAMPLE" << sampleIndex << " " << li.getSample(sampleIndex);
    }
    return os;
}



std::ostream&
operator<<(
    std::ostream& os,
    const GermlineIndelLocusInfo& ii)
{
    ii.assertValidity();

    os << "GermlineIndelLocusInfo ";
    os << static_cast<const LocusInfo&>(ii);

    const auto& altAlleles(ii.getIndelAlleles());
    os << "AltAlleleCount: " << altAlleles.size() << "\n";
    for (unsigned altAlleleIndex(0); altAlleleIndex<altAlleles.size(); ++altAlleleIndex)
    {
        os << "AltAllele index/info: " << altAlleleIndex << " " << altAlleles[altAlleleIndex];
    }
    os << "range: " << ii.range() << "\n";
    const unsigned sampleCount(ii.getSampleCount());
    for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
    {
        const auto& indelSampleInfo(ii.getIndelSample(sampleIndex));
        os << "IndelSample" << sampleIndex << "\n";
        os << indelSampleInfo.legacyReportInfo << "\n";
    }

    return os;
}



std::ostream&
operator<<(
    std::ostream& os,
    const GermlineSiteLocusInfo& si)
{
    os << "GermlineSiteLocusInfo ";
    os << static_cast<const LocusInfo&>(si);

    os << "ref: " << id_to_base(si.refBaseIndex) << "\n";

    const auto& altAlleles(si.getSiteAlleles());
    os << "AltAlleleCount: " << altAlleles.size() << "\n";
    for (unsigned altAlleleIndex(0); altAlleleIndex<altAlleles.size(); ++altAlleleIndex)
    {
        os << "AltAllele index/info: " << altAlleleIndex << " " << altAlleles[altAlleleIndex] << "\n";
    }

    return os;
}



std::ostream&
operator<<(
    std::ostream& os,
    const GermlineDiploidSiteLocusInfo& si)
{
    os << "GermlineDiploidSiteLocusInfo ";
    os << static_cast<const GermlineSiteLocusInfo&>(si);

    return os;
}



/// code variant as expected in the EVS model
static
unsigned
getEVSGenotypeCode(
    const bool isDiploid,
    const uint8_t allele0Index,
    const uint8_t allele1Index)
{
    enum {HET, HOM, HETALT};

    if (isDiploid)
    {
        if (allele0Index == allele1Index)
        {
            return HOM;
        }
        else
        {
            if ((allele0Index > 0) and (allele1Index > 0))
            {
                return HETALT;
            }
            else
            {
                return HET;
            }
        }
    }
    else
    {
        return HOM;
    }
}



static
void
getAlleleBiasFeatures(
    const LocusSupportingReadStats& counts,
    const bool isDiploid,
    const uint8_t allele0Index,
    const uint8_t allele1Index,
    double& SampleIndelAlleleBiasLower,
    double& SampleIndelAlleleBias)
{
    // alleleBiasLower cdf of binomial prob of seeing no more than the number of 'allele A' reads out of A reads + B reads, given p=0.5
    // alleleBiasUpper cdf of binomial prob of seeing no more than the number of 'allele B' reads out of A reads + B reads, given p=0.5
    double alleleBiasLower(0);
    double alleleBiasUpper(0);

    if (isDiploid and (allele0Index != allele1Index))
    {
        const auto& fwdCounts(counts.getCounts(true));
        const auto& revCounts(counts.getCounts(false));

        const double r0(fwdCounts.confidentAlleleCount(allele0Index) + revCounts.confidentAlleleCount(allele0Index));
        const double r1(fwdCounts.confidentAlleleCount(allele1Index) + revCounts.confidentAlleleCount(allele1Index));
        alleleBiasLower = cdf(boost::math::binomial(r0 + r1, 0.5), r0);
        alleleBiasUpper = cdf(boost::math::binomial(r0 + r1, 0.5), r1);
    }

    // fudge to avoid log(0) in extreme cases
    static const double fudge(1.e-30);
    SampleIndelAlleleBiasLower = -std::log(alleleBiasLower + fudge);
    SampleIndelAlleleBias =-std::log(std::min(1., 2. * std::min(alleleBiasLower, alleleBiasUpper)) + fudge);
}


void
GermlineDiploidSiteLocusInfo::
computeEmpiricalScoringFeatures(
    const GermlineDiploidSiteLocusInfo& locus,
    const unsigned sampleIndex,
    const bool isRNA,
    const bool isUniformDepthExpected,
    const bool isComputeDevelopmentFeatures,
    const double allSampleChromDepth,
    VariantScoringFeatureKeeper& features,
    VariantScoringFeatureKeeper& developmentFeatures)
{
    const auto& sampleInfo(locus.getSample(sampleIndex));
    const auto& siteSampleInfo(locus.getSiteSample(sampleIndex));

    // filtered locus depth/total locus depth/confidentDepth for this sample:
    const double filteredLocusDepth(siteSampleInfo.usedBasecallCount);
    const double totalLocusDepth(siteSampleInfo.getTotalReadDepth());

    // total locus depth summed over all samples and MAPQ values
    // (this makes it comparable to chromDepth, which is computed over all samples/MAPQ):
    /// TODO account for local copy number in both total locus depth and expected (chrom) depth
    const double allSampleTotalLocusDepth(locus.getTotalReadDepth());

    const double allSampleChromDepthFactor(safeFrac(1, allSampleChromDepth));
    const double filteredLocusDepthFactor(safeFrac(1, filteredLocusDepth));
    const double totalLocusDepthFactor(safeFrac(1, totalLocusDepth));

    const auto& siteAlleles(locus.getSiteAlleles());
    assert(not siteAlleles.empty());

    uint8_t allele0Index(0);
    uint8_t allele1Index(0);
    unsigned primaryAltAlleleIndex(0);

    const auto& maxGt(sampleInfo.max_gt());
    const SamplePloidyState ploidy(maxGt.getPloidy());

    if (ploidy.isDiploid())
    {
        allele0Index = maxGt.getAllele0Index();
        allele1Index = maxGt.getAllele1Index();
        if (allele0Index>0)
        {
            primaryAltAlleleIndex = (allele0Index - 1);
        }
        else if (allele1Index>0)
        {
            primaryAltAlleleIndex = (allele1Index - 1);
        }
    }
    else if (ploidy.isHaploid())
    {
        allele0Index = maxGt.getAllele0Index();
        if (allele0Index>0)
        {
            primaryAltAlleleIndex = allele0Index-1;
        }
    }
    else
    {
        assert(false and "Unexpected ploidy");
    }

    const unsigned confidentRefCount(
        sampleInfo.supportCounts.getCounts(true).confidentRefAlleleCount() +
        sampleInfo.supportCounts.getCounts(false).confidentRefAlleleCount());
    const unsigned confidentPrimaryAltCount(
        sampleInfo.supportCounts.getCounts(true).confidentAltAlleleCount(primaryAltAlleleIndex) +
        sampleInfo.supportCounts.getCounts(false).confidentAltAlleleCount(primaryAltAlleleIndex));

    const double mapqZeroFraction(siteSampleInfo.mapqTracker.getZeroFrac());

    const double locusUsedDepthFraction(filteredLocusDepth * totalLocusDepthFactor);


    if (isRNA)
    {
        features.set(RNA_SNV_SCORING_FEATURES::SiteHomopolymerLength, (locus.hpol));
        features.set(RNA_SNV_SCORING_FEATURES::SampleStrandBias, (siteSampleInfo.strandBias));
        features.set(RNA_SNV_SCORING_FEATURES::SamplePrimaryAltAlleleDepth, (confidentPrimaryAltCount));
        features.set(RNA_SNV_SCORING_FEATURES::SamplePrimaryAltAlleleDepthFraction, (confidentPrimaryAltCount * filteredLocusDepthFactor));
        features.set(RNA_SNV_SCORING_FEATURES::VariantAlleleQuality, (locus.anyVariantAlleleQuality));
        features.set(RNA_SNV_SCORING_FEATURES::SampleMeanDistanceFromReadEdge, (siteSampleInfo.meanDistanceFromReadEdge));

        // compute any experimental features not currently used in production
        //
        if (isComputeDevelopmentFeatures)
        {
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::GT, getEVSGenotypeCode(ploidy.isDiploid(), allele0Index, allele1Index));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::SampleRefAlleleDepth, (confidentRefCount));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::F_DP, (siteSampleInfo.usedBasecallCount));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::F_DPF, (siteSampleInfo.unusedBasecallCount));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::I_BaseQRankSum, (siteSampleInfo.BaseQRankSum));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::SampleReadPosRankSum, (siteSampleInfo.ReadPosRankSum));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::I_AvgBaseQ, (siteSampleInfo.avgBaseQ));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::SampleRMSMappingQuality, (siteSampleInfo.mapqTracker.getRMS()));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::SampleRMSMappingQualityRankSum, (siteSampleInfo.MQRankSum));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::mapqZeroFraction, mapqZeroFraction);
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::SampleUsedDepthFraction, locusUsedDepthFraction);
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::QUAL_NORM,
                                    (locus.anyVariantAlleleQuality * filteredLocusDepthFactor));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQX_NORM,
                                    (sampleInfo.gqx * filteredLocusDepthFactor));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQ_NORM,
                                    (sampleInfo.genotypeQualityPolymorphic * filteredLocusDepthFactor));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::SampleRefAlleleDepthFraction,
                                    (confidentRefCount * filteredLocusDepthFactor));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::ConservativeGenotypeQuality, (sampleInfo.gqx));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQ, (sampleInfo.genotypeQualityPolymorphic));
        }
    }
    else
    {
        features.set(GERMLINE_SNV_SCORING_FEATURES::GenotypeCategory,
                     getEVSGenotypeCode(ploidy.isDiploid(), allele0Index, allele1Index));

        features.set(GERMLINE_SNV_SCORING_FEATURES::SampleRMSMappingQuality, (siteSampleInfo.mapqTracker.getRMS()));
        features.set(GERMLINE_SNV_SCORING_FEATURES::SiteHomopolymerLength, (locus.hpol));
        features.set(GERMLINE_SNV_SCORING_FEATURES::SampleStrandBias, (siteSampleInfo.strandBias));
        features.set(GERMLINE_SNV_SCORING_FEATURES::SampleRMSMappingQualityRankSum, (siteSampleInfo.MQRankSum));
        features.set(GERMLINE_SNV_SCORING_FEATURES::SampleReadPosRankSum, (siteSampleInfo.ReadPosRankSum));

        // how surprising is the depth relative to expect? This is the only EVS feature modified for exome/targeted runs
        /// TODO: convert this to pvalue based on Poisson distro?
        double relativeLocusDepth(1.);
        if (isUniformDepthExpected)
        {
            relativeLocusDepth = (allSampleTotalLocusDepth * allSampleChromDepthFactor);
        }

        features.set(GERMLINE_SNV_SCORING_FEATURES::RelativeTotalLocusDepth, relativeLocusDepth);

        // how noisy is the locus?
        features.set(GERMLINE_SNV_SCORING_FEATURES::SampleUsedDepthFraction, locusUsedDepthFraction);

        features.set(GERMLINE_SNV_SCORING_FEATURES::ConservativeGenotypeQuality, (sampleInfo.gqx));

        // how confident are we in the haplotyping step?
        float normalizedAltHaplotypeCountRatio;
        if (locus.getActiveRegionId() < 0)
            normalizedAltHaplotypeCountRatio = -1.0f;   // not in active region
        else if (maxGt.getPloidy() == 1)
            normalizedAltHaplotypeCountRatio = maxGt.getAltHaplotypeCountRatio()*2.0f; // multiply by 2 for ploidy=1
        else if (maxGt.isHet() and maxGt.containsReference())
            normalizedAltHaplotypeCountRatio = maxGt.getAltHaplotypeCountRatio()*2.0f; // multiply by 2 for het (but not hetalt)
        else
            normalizedAltHaplotypeCountRatio = maxGt.getAltHaplotypeCountRatio();
        features.set(GERMLINE_SNV_SCORING_FEATURES::NormalizedAltHaplotypeCountRatio, normalizedAltHaplotypeCountRatio);

        // compute any experimental features not currently used in production
        //
        if (isComputeDevelopmentFeatures)
        {
            // BaseQRankSum
            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::I_BaseQRankSum, (siteSampleInfo.BaseQRankSum));

            // allele bias metrics
            {
                double SampleIndelAlleleBiasLower, SampleIndelAlleleBias;
                getAlleleBiasFeatures(
                    sampleInfo.supportCounts, ploidy.isDiploid(), allele0Index, allele1Index,
                    SampleIndelAlleleBiasLower, SampleIndelAlleleBias);

                developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::SampleIndelAlleleBiasLower, SampleIndelAlleleBiasLower);
                developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::SampleIndelAlleleBias, SampleIndelAlleleBias);
            }

            //The average baseQ of the position of alt allele
            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::I_RawBaseQ, (siteSampleInfo.avgBaseQ));

            //the average position value within a read of alt allele
            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::I_RawPos, (siteSampleInfo.rawPos));

            // hom unreliable are the read mappings near this locus?
            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::mapqZeroFraction, mapqZeroFraction);

            // renormalized features intended to replace the corresponding production feature
            //
            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::QUAL_NORM,
                                    (locus.anyVariantAlleleQuality * filteredLocusDepthFactor));
            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQX_NORM,
                                    (sampleInfo.gqx * filteredLocusDepthFactor));
            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQ_NORM,
                                    (sampleInfo.genotypeQualityPolymorphic * filteredLocusDepthFactor));

            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::SampleRefAlleleDepthFraction,
                                    (confidentRefCount * filteredLocusDepthFactor));

            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::VariantAlleleQuality,
                                    (locus.anyVariantAlleleQuality));
            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQ, (sampleInfo.genotypeQualityPolymorphic));

            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::SamplePrimaryAltAlleleDepthFraction,
                                    (confidentPrimaryAltCount * filteredLocusDepthFactor));
        }
    }
}



void
GermlineDiploidIndelLocusInfo::
computeEmpiricalScoringFeatures(
    const GermlineDiploidIndelLocusInfo& locus,
    const unsigned sampleIndex,
    const bool isRNA,
    const bool isUniformDepthExpected,
    const bool isComputeDevelopmentFeatures,
    const double allSampleChromDepth,
    VariantScoringFeatureKeeper& features,
    VariantScoringFeatureKeeper& developmentFeatures)
{
    const auto& sampleInfo(locus.getSample(sampleIndex));
    const auto& indelSampleInfo(locus.getIndelSample(sampleIndex));

    // filtered locus depth/total locus depth/confidentDepth for this sample:
    const double filteredLocusDepth(indelSampleInfo.tier1Depth);
    const double locusDepth(indelSampleInfo.getTotalReadDepth());
    const double confidentDepth(sampleInfo.supportCounts.totalConfidentCounts());

    // total locus depth summed over all samples (compatible with chromDepth, which is computed over all samples):
    /// TODO account for local copy number in both total locus depth and expected (chrom) depth
    const double allSampleLocusDepth(locus.getTotalReadDepth());

    const double allSampleChromDepthFactor(safeFrac(1,allSampleChromDepth));
    const double filteredLocusDepthFactor(safeFrac(1,filteredLocusDepth));
    const double locusDepthFactor(safeFrac(1,locusDepth));
    const double confidentDepthFactor(safeFrac(1,confidentDepth));

    uint8_t allele0Index(0);
    uint8_t allele1Index(0);
    unsigned primaryAltAlleleIndex(0);

    const auto& maxGt(sampleInfo.max_gt());
    const SamplePloidyState ploidy(maxGt.getPloidy());

    if (ploidy.isDiploid())
    {
        allele0Index = maxGt.getAllele0Index();
        allele1Index = maxGt.getAllele1Index();
        if (allele0Index>0)
        {
            primaryAltAlleleIndex = (allele0Index - 1);
        }
        else if (allele1Index>0)
        {
            primaryAltAlleleIndex = (allele1Index - 1);
        }
    }
    else if (ploidy.isHaploid())
    {
        allele0Index = maxGt.getAllele0Index();
        if (allele0Index>0)
        {
            primaryAltAlleleIndex = allele0Index-1;
        }
    }
    else
    {
        assert(false and "Unexpected ploidy");
    }

    ///TODO STREL-125 generalize to multiple alts:
    const auto& primaryAltAllele(locus.getIndelAlleles()[primaryAltAlleleIndex]);

    const unsigned confidentRefCount(
        sampleInfo.supportCounts.getCounts(true).confidentRefAlleleCount() +
        sampleInfo.supportCounts.getCounts(false).confidentRefAlleleCount());
    const unsigned confidentPrimaryAltFwdCount(
        sampleInfo.supportCounts.fwdCounts.confidentAltAlleleCount(primaryAltAlleleIndex));
    const unsigned confidentPrimaryAltRevCount(
        sampleInfo.supportCounts.revCounts.confidentAltAlleleCount(primaryAltAlleleIndex));
    const unsigned confidentPrimaryAltCount(confidentPrimaryAltFwdCount + confidentPrimaryAltRevCount);

    // strand bias from allele counts
    const unsigned totalConfidentFwdCount(sampleInfo.supportCounts.fwdCounts.totalConfidentCounts());
    const unsigned totalConfidentRevCount(sampleInfo.supportCounts.revCounts.totalConfidentCounts());
    const double alleleCountRawStrandBias(starling_continuous_variant_caller::strandBias(
                                              confidentPrimaryAltFwdCount,
                                              confidentPrimaryAltRevCount,
                                              totalConfidentFwdCount - confidentPrimaryAltFwdCount,
                                              totalConfidentRevCount - confidentPrimaryAltRevCount));

    // Bound 'raw' strand-bias input to be less than the specified absolute value
    static const double EVSMaxAbsSampleVariantStrandBias(99.);
    const double alleleCountStrandBias(std::min(EVSMaxAbsSampleVariantStrandBias, std::max(-EVSMaxAbsSampleVariantStrandBias, alleleCountRawStrandBias)));

    // allele bias metrics
    double SampleIndelAlleleBiasLower, SampleIndelAlleleBias;
    getAlleleBiasFeatures(sampleInfo.supportCounts, ploidy.isDiploid(), allele0Index, allele1Index, SampleIndelAlleleBiasLower, SampleIndelAlleleBias);

    if (isRNA)
    {
        features.set(RNA_INDEL_SCORING_FEATURES::SampleRefRepeatCount, (primaryAltAllele.indelReportInfo.refRepeatCount));
        features.set(RNA_INDEL_SCORING_FEATURES::SampleIndelRepeatCount, (primaryAltAllele.indelReportInfo.indelRepeatCount));
        features.set(RNA_INDEL_SCORING_FEATURES::SampleIndelRepeatUnitSize, (primaryAltAllele.indelReportInfo.repeatUnit.length()));
        features.set(RNA_INDEL_SCORING_FEATURES::SampleRefAlleleDepth, (confidentRefCount));
        features.set(RNA_INDEL_SCORING_FEATURES::SamplePrimaryAltAlleleDepth, (confidentPrimaryAltCount));
        features.set(RNA_INDEL_SCORING_FEATURES::SamplePrimaryAltAlleleDepthFraction, (confidentPrimaryAltCount * confidentDepthFactor));
        features.set(RNA_INDEL_SCORING_FEATURES::VariantAlleleQuality, (locus.anyVariantAlleleQuality));
        // note this feature's behavior hasn't been figured out for an overlapping indel locus
        features.set(RNA_INDEL_SCORING_FEATURES::SampleIndelMeanDistanceFromReadEdge, (indelSampleInfo.legacyReportInfo.distanceFromReadEdge.mean()));

        // compute any experimental features not currently used in production
        if (isComputeDevelopmentFeatures)
        {
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_DPI, (indelSampleInfo.tier1Depth));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::SampleIndelAlleleBiasLower, SampleIndelAlleleBiasLower);
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::SampleIndelAlleleBias, SampleIndelAlleleBias);

            // how unreliable are the read mappings near this locus?
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::SampleProxyRMSMappingQuality,
                                    (indelSampleInfo.mapqTracker.getRMS()));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::mapqZeroFraction,
                                    (indelSampleInfo.mapqTracker.getZeroFrac()));

            // how noisy is the locus?
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_DPI_NORM,
                                    (filteredLocusDepth * locusDepthFactor));

            // all of the features below are simply renormalized replacements of the current production feature set
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::QUAL_NORM,
                                    (locus.anyVariantAlleleQuality * filteredLocusDepthFactor));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQX_NORM,
                                    (sampleInfo.gqx * filteredLocusDepthFactor));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQ_NORM,
                                    (sampleInfo.genotypeQualityPolymorphic * filteredLocusDepthFactor));

            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::SampleRefAlleleDepthFraction,
                                    (confidentRefCount * confidentDepthFactor));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::ConservativeGenotypeQuality, (sampleInfo.gqx));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQ, (sampleInfo.genotypeQualityPolymorphic));
        }
    }
    else
    {
        features.set(GERMLINE_INDEL_SCORING_FEATURES::GenotypeCategory,
                     getEVSGenotypeCode(ploidy.isDiploid(), allele0Index, allele1Index));
        features.set(GERMLINE_INDEL_SCORING_FEATURES::SampleIndelRepeatCount, (primaryAltAllele.indelReportInfo.indelRepeatCount));
        features.set(GERMLINE_INDEL_SCORING_FEATURES::SampleIndelRepeatUnitSize, (primaryAltAllele.indelReportInfo.repeatUnit.length()));

        features.set(GERMLINE_INDEL_SCORING_FEATURES::SampleIndelAlleleBiasLower, SampleIndelAlleleBiasLower);
        features.set(GERMLINE_INDEL_SCORING_FEATURES::SampleIndelAlleleBias, SampleIndelAlleleBias);

        // how unreliable are the read mappings near this locus?
        features.set(GERMLINE_INDEL_SCORING_FEATURES::SampleProxyRMSMappingQuality,
                     (indelSampleInfo.mapqTracker.getRMS()));

        // how surprising is the depth relative to expect? This is the only EVS feature modified for exome/targeted runs
        //
        /// TODO: convert this to pvalue based on Poisson distro?
        double relativeLocusDepth(1.);
        if (isUniformDepthExpected)
        {
            relativeLocusDepth = (allSampleLocusDepth * allSampleChromDepthFactor);
        }
        features.set(GERMLINE_INDEL_SCORING_FEATURES::RelativeTotalLocusDepth, relativeLocusDepth);

        features.set(GERMLINE_INDEL_SCORING_FEATURES::SamplePrimaryAltAlleleDepthFraction,
                     (confidentPrimaryAltCount * confidentDepthFactor));

        features.set(GERMLINE_INDEL_SCORING_FEATURES::ConservativeGenotypeQuality, (sampleInfo.gqx));
        features.set(GERMLINE_INDEL_SCORING_FEATURES::InterruptedHomopolymerLength, (primaryAltAllele.indelReportInfo.interruptedHomopolymerLength));
        features.set(GERMLINE_INDEL_SCORING_FEATURES::ContextCompressability, (primaryAltAllele.indelReportInfo.contextCompressability));
        features.set(GERMLINE_INDEL_SCORING_FEATURES::IndelCategory, (primaryAltAllele.indelKey.isPrimitiveDeletionAllele()));
        features.set(GERMLINE_INDEL_SCORING_FEATURES::SampleAlleleCountStrandBias, (alleleCountStrandBias));

        // how confident are we in the haplotyping step?
        float normalizedAltHaplotypeCountRatio;
        if (locus.getActiveRegionId() < 0)
            normalizedAltHaplotypeCountRatio = -1.0f;   // not in active region
        else if (maxGt.getPloidy() == 1)
            normalizedAltHaplotypeCountRatio = maxGt.getAltHaplotypeCountRatio()*2.0f; // multiply by 2 for ploidy=1
        else if (maxGt.isHet() and maxGt.containsReference())
            normalizedAltHaplotypeCountRatio = maxGt.getAltHaplotypeCountRatio()*2.0f; // multiply by 2 for het (but not hetalt)
        else
            normalizedAltHaplotypeCountRatio = maxGt.getAltHaplotypeCountRatio();
        features.set(GERMLINE_INDEL_SCORING_FEATURES::NormalizedAltHaplotypeCountRatio, normalizedAltHaplotypeCountRatio);

        // compute any experimental features not currently used in production
        if (isComputeDevelopmentFeatures)
        {
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::SampleRefRepeatCount, (primaryAltAllele.indelReportInfo.refRepeatCount));

            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::mapqZeroFraction,
                                    (indelSampleInfo.mapqTracker.getZeroFrac()));

            // how noisy is the locus?
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::F_DPI_NORM,
                                    (filteredLocusDepth * locusDepthFactor));

            // all of the features below are simply renormalized replacements of the current production feature set
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::QUAL_NORM,
                                    (locus.anyVariantAlleleQuality * filteredLocusDepthFactor));
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQX_NORM,
                                    (sampleInfo.gqx * filteredLocusDepthFactor));
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQ_NORM,
                                    (sampleInfo.genotypeQualityPolymorphic * filteredLocusDepthFactor));

            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::SampleRefAlleleDepthFraction,
                                    (confidentRefCount * confidentDepthFactor));

            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::VariantAlleleQuality, (locus.anyVariantAlleleQuality));
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQ, (sampleInfo.genotypeQualityPolymorphic));
        }
    }
}
