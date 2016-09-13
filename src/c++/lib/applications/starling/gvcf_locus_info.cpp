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

///
/// \author Chris Saunders
///


#include "gvcf_locus_info.hh"
#include "blt_util/math_util.hh"
#include "common/Exceptions.hh"

#include "boost/math/distributions/binomial.hpp"
#include "boost/math/distributions.hpp"
#include "rnaVariantEmpiricalScoringFeatures.hh"

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



void
GermlineIndelLocusInfo::
getOffsetError(
    const unsigned offset) const
{
    using namespace illumina::common;

    std::ostringstream oss;
    oss << "ERROR: indel locus site offset '" << offset << "' exceeds exceeds locus size '" << _range.size() << "'\n";
    BOOST_THROW_EXCEPTION(LogicException(oss.str()));
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
operator<<(std::ostream& os,
           const GermlineSiteLocusInfo& si)
{
    os << "GermlineSiteLocusInfo ";
    os << static_cast<const LocusInfo&>(si);
    return os;
}



std::ostream&
operator<<(std::ostream& os,
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
    double& ABlower,
    double& AB)
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
    ABlower = -std::log(alleleBiasLower + fudge);
    AB =-std::log(std::min(1., 2. * std::min(alleleBiasLower, alleleBiasUpper)) + fudge);
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
    const double filteredLocusDepth(siteSampleInfo.n_used_calls);
    const double locusDepth(siteSampleInfo.getTotalReadDepth());

    // total locus depth summed over all samples (compatible with chromDepth, which is computed over all samples):
    /// TODO account for local copy number in both total locus depth and expected (chrom) depth
    const double allSampleLocusDepth(locus.getTotalReadDepth());

    const double chromDepthFactor(safeFrac(1, allSampleChromDepth));
    const double filteredLocusDepthFactor(safeFrac(1, filteredLocusDepth));
    const double locusDepthFactor(safeFrac(1, locusDepth));

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

    ///TODO STREL-125 generalize to multiple alts:
    const auto& primaryAltAllele(siteAlleles[primaryAltAlleleIndex]);

    const unsigned confidentRefCount(
        sampleInfo.supportCounts.getCounts(true).confidentRefAlleleCount() +
        sampleInfo.supportCounts.getCounts(false).confidentRefAlleleCount());
    const unsigned confidentPrimaryAltCount(
        sampleInfo.supportCounts.getCounts(true).confidentAltAlleleCount(primaryAltAlleleIndex) +
        sampleInfo.supportCounts.getCounts(false).confidentAltAlleleCount(primaryAltAlleleIndex));

    const double mapqZeroFraction(siteSampleInfo.mapqTracker.getZeroFrac());

    const double locusUsedDepthFraction(filteredLocusDepth * locusDepthFactor);


    if (isRNA)
    {
        features.set(RNA_SNV_SCORING_FEATURES::GT,
                     getEVSGenotypeCode(ploidy.isDiploid(), allele0Index, allele1Index));

        features.set(RNA_SNV_SCORING_FEATURES::QUAL, (locus.anyVariantAlleleQuality * chromDepthFactor));
        features.set(RNA_SNV_SCORING_FEATURES::F_DP, (siteSampleInfo.n_used_calls * chromDepthFactor));
        features.set(RNA_SNV_SCORING_FEATURES::F_DPF, (siteSampleInfo.n_unused_calls * chromDepthFactor));
        features.set(RNA_SNV_SCORING_FEATURES::F_GQ, (sampleInfo.genotypeQualityPolymorphic * chromDepthFactor));
        features.set(RNA_SNV_SCORING_FEATURES::F_GQX, (sampleInfo.gqx * chromDepthFactor));

        features.set(RNA_SNV_SCORING_FEATURES::I_AvgBaseQ, (siteSampleInfo.avgBaseQ));
        features.set(RNA_SNV_SCORING_FEATURES::I_AvgPos, (siteSampleInfo.rawPos));

        features.set(RNA_SNV_SCORING_FEATURES::I_BaseQRankSum, (siteSampleInfo.BaseQRankSum));
        features.set(RNA_SNV_SCORING_FEATURES::I_ReadPosRankSum, (siteSampleInfo.ReadPosRankSum));

        features.set(RNA_SNV_SCORING_FEATURES::I_SNVHPOL, (locus.hpol));
        features.set(RNA_SNV_SCORING_FEATURES::I_SNVSB, (primaryAltAllele.strandBias));

        features.set(RNA_SNV_SCORING_FEATURES::AD0, (confidentRefCount * chromDepthFactor));
        features.set(RNA_SNV_SCORING_FEATURES::AD1, (confidentPrimaryAltCount * chromDepthFactor));

        /// TODO STREL-125 generalize this to multi-alts:
        features.set(RNA_SNV_SCORING_FEATURES::ADR, safeFrac(confidentRefCount, (confidentRefCount + confidentPrimaryAltCount)));

        // compute any experimental features not currently used in production
        //
        if (isComputeDevelopmentFeatures)
        {
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::I_MQ, (siteSampleInfo.mapqTracker.getRMS()));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::I_MQRankSum, (siteSampleInfo.MQRankSum));

            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::mapqZeroFraction, mapqZeroFraction);

            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::F_DP_NORM, locusUsedDepthFraction);

            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::QUAL_NORM,
                                          (locus.anyVariantAlleleQuality * filteredLocusDepthFactor));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQX_NORM,
                                          (sampleInfo.gqx * filteredLocusDepthFactor));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQ_NORM,
                                          (sampleInfo.genotypeQualityPolymorphic * filteredLocusDepthFactor));

            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::AD0_NORM,
                                          (confidentRefCount * filteredLocusDepthFactor));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::AD1_NORM,
                                          (confidentPrimaryAltCount * filteredLocusDepthFactor));

            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::QUAL_EXACT,
                                          (locus.anyVariantAlleleQuality));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQX_EXACT, (sampleInfo.gqx));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQ_EXACT, (sampleInfo.genotypeQualityPolymorphic));
        }
    }
    else
    {
        features.set(GERMLINE_SNV_SCORING_FEATURES::GENO,
                     getEVSGenotypeCode(ploidy.isDiploid(), allele0Index, allele1Index));

        features.set(GERMLINE_SNV_SCORING_FEATURES::I_MQ, (siteSampleInfo.mapqTracker.getRMS()));
        features.set(GERMLINE_SNV_SCORING_FEATURES::I_SNVHPOL, (locus.hpol));
        features.set(GERMLINE_SNV_SCORING_FEATURES::I_SNVSB, (primaryAltAllele.strandBias));
        features.set(GERMLINE_SNV_SCORING_FEATURES::I_MQRankSum, (siteSampleInfo.MQRankSum));
        features.set(GERMLINE_SNV_SCORING_FEATURES::I_ReadPosRankSum, (siteSampleInfo.ReadPosRankSum));

        // how surprising is the depth relative to expect? This is the only value will be modified for exome/targeted runs
        /// TODO: convert this to pvalue based on Poisson distro?
        double relativeLocusDepth(1.);
        if (isUniformDepthExpected)
        {
            relativeLocusDepth = (allSampleLocusDepth * chromDepthFactor);
        }

        features.set(GERMLINE_SNV_SCORING_FEATURES::TDP_NORM, relativeLocusDepth);

        // how noisy is the locus?
        features.set(GERMLINE_SNV_SCORING_FEATURES::F_DP_NORM, locusUsedDepthFraction);

        features.set(GERMLINE_SNV_SCORING_FEATURES::F_GQX_EXACT, (sampleInfo.gqx));

        // compute any experimental features not currently used in production
        //
        if (isComputeDevelopmentFeatures)
        {

            // BaseQRankSum
            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::I_BaseQRankSum, (siteSampleInfo.BaseQRankSum));

            // allele bias metrics
            {
                double ABlower, AB;
                getAlleleBiasFeatures(
                    sampleInfo.supportCounts, ploidy.isDiploid(), allele0Index, allele1Index,
                    ABlower, AB);

                developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::ABlower, ABlower);
                developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::AB, AB);
            }

            //The average baseQ of the position of alt allele
            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::I_RawBaseQ, (siteSampleInfo.avgBaseQ));

            //the average position value within a read of alt allele
            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::I_RawPos, (siteSampleInfo.rawPos));

            // hom unrelable are the read mappings near this locus?
            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::mapqZeroFraction, mapqZeroFraction);

            // renormalized features intended to replace the corresponding production feature
            //
            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::QUAL_NORM,
                                          (locus.anyVariantAlleleQuality * filteredLocusDepthFactor));
            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQX_NORM,
                                          (sampleInfo.gqx * filteredLocusDepthFactor));
            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQ_NORM,
                                          (sampleInfo.genotypeQualityPolymorphic * filteredLocusDepthFactor));

            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::AD0_NORM,
                                          (confidentRefCount * filteredLocusDepthFactor));

            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::QUAL_EXACT,
                                          (locus.anyVariantAlleleQuality));
            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQ_EXACT, (sampleInfo.genotypeQualityPolymorphic));

            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::AD1_NORM,
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
    const unsigned confidentPrimaryAltCount(
        sampleInfo.supportCounts.getCounts(true).confidentAltAlleleCount(primaryAltAlleleIndex) +
        sampleInfo.supportCounts.getCounts(false).confidentAltAlleleCount(primaryAltAlleleIndex));

    // allele bias metrics
    double ABlower, AB;
    getAlleleBiasFeatures(sampleInfo.supportCounts, ploidy.isDiploid(), allele0Index, allele1Index, ABlower, AB);

        if (isRNA)
    {
        features.set(RNA_INDEL_SCORING_FEATURES::QUAL, (locus.anyVariantAlleleQuality * allSampleChromDepthFactor));
        features.set(RNA_INDEL_SCORING_FEATURES::F_GQX, (sampleInfo.gqx * allSampleChromDepthFactor));
        features.set(RNA_INDEL_SCORING_FEATURES::REFREP1, (primaryAltAllele.indelReportInfo.ref_repeat_count));
        features.set(RNA_INDEL_SCORING_FEATURES::IDREP1, (primaryAltAllele.indelReportInfo.indel_repeat_count));
        features.set(RNA_INDEL_SCORING_FEATURES::RULEN1, (primaryAltAllele.indelReportInfo.repeat_unit.length()));
        features.set(RNA_INDEL_SCORING_FEATURES::AD0,
                     (confidentRefCount * allSampleChromDepthFactor));
        features.set(RNA_INDEL_SCORING_FEATURES::AD1,
                     (confidentPrimaryAltCount * allSampleChromDepthFactor));
        features.set(RNA_INDEL_SCORING_FEATURES::F_DPI, (indelSampleInfo.tier1Depth * allSampleChromDepthFactor));

        developmentFeatures.set(RNA_INDEL_SCORING_FEATURES::ABlower, ABlower);
        developmentFeatures.set(RNA_INDEL_SCORING_FEATURES::AB, AB);

        // compute any experimental features not currently used in production
        if (isComputeDevelopmentFeatures)
        {
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQ, (sampleInfo.genotypeQualityPolymorphic * allSampleChromDepthFactor));

            // how unreliable are the read mappings near this locus?
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_MQ,
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

            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::AD0_NORM,
                                    (confidentRefCount * confidentDepthFactor));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::AD1_NORM,
                                    (confidentPrimaryAltCount * confidentDepthFactor));

            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::QUAL_EXACT, (locus.anyVariantAlleleQuality));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQX_EXACT, (sampleInfo.gqx));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQ_EXACT, (sampleInfo.genotypeQualityPolymorphic));
        }
    }
    else
    {
        features.set(GERMLINE_INDEL_SCORING_FEATURES::GENO,
                     getEVSGenotypeCode(ploidy.isDiploid(), allele0Index, allele1Index));
        features.set(GERMLINE_INDEL_SCORING_FEATURES::IDREP1, (primaryAltAllele.indelReportInfo.indel_repeat_count));
        features.set(GERMLINE_INDEL_SCORING_FEATURES::RULEN1, (primaryAltAllele.indelReportInfo.repeat_unit.length()));

        features.set(GERMLINE_INDEL_SCORING_FEATURES::ABlower, ABlower);
        features.set(GERMLINE_INDEL_SCORING_FEATURES::AB, AB);

        // how unreliable are the read mappings near this locus?
        features.set(GERMLINE_INDEL_SCORING_FEATURES::F_MQ,
                     (indelSampleInfo.mapqTracker.getRMS()));

        // how surprising is the depth relative to expect? This is the only value will be modified for exome/targeted runs
        //
        /// TODO: convert this to pvalue based on Poisson distro?
        double relativeLocusDepth(1.);
        if (isUniformDepthExpected)
        {
            relativeLocusDepth = (allSampleLocusDepth * allSampleChromDepthFactor);
        }
        features.set(GERMLINE_INDEL_SCORING_FEATURES::TDP_NORM, relativeLocusDepth);

        features.set(GERMLINE_INDEL_SCORING_FEATURES::AD1_NORM,
                     (confidentPrimaryAltCount * confidentDepthFactor));

        features.set(GERMLINE_INDEL_SCORING_FEATURES::F_GQX_EXACT, (sampleInfo.gqx));

        // compute any experimental features not currently used in production
        if (isComputeDevelopmentFeatures)
        {
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::REFREP1, (primaryAltAllele.indelReportInfo.ref_repeat_count));

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

            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::AD0_NORM,
                                    (confidentRefCount * confidentDepthFactor));

            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::QUAL_EXACT, (locus.anyVariantAlleleQuality));
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQ_EXACT, (sampleInfo.genotypeQualityPolymorphic));
        }
    }
}
