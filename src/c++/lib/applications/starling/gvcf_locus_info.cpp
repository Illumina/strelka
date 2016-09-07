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
       << " gqs: " << lsi.gqx
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
        os << "AltAllele" << altAlleleIndex << " " << altAlleles[altAlleleIndex].indelKey;
    }
    os << "range: " << ii.range() << "\n";
    const unsigned sampleCount(ii.getSampleCount());
    for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
    {
        const auto& indelSampleInfo(ii.getIndelSample(sampleIndex));
        os << "IndelSample" << sampleIndex << "\n";
        os << indelSampleInfo.legacyReportInfo << "\n";
        os << "sp size: " << indelSampleInfo.sitePloidy.size() << "\n";
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

    // get the alt base id (choose second in case of an alt het....)
    unsigned altBase(N_BASE);
    const uint8_t refBaseId(base_to_id(locus.ref));
    for (unsigned b(0); b < N_BASE; ++b)
    {
        if (b == refBaseId) continue;
        if (DIGT::expect2(b, siteSampleInfo.max_gt))
        {
            altBase = b;
        }
    }
    assert(altBase != N_BASE);

    /// TODO STREL-125 generalize to multi-alt
    const unsigned r0 = siteSampleInfo.alleleObservationCounts(refBaseId);
    const unsigned r1 = siteSampleInfo.alleleObservationCounts(altBase);

    const double mapqZeroFraction(safeFrac(locus.mapqZeroCount, locus.mapqCount));

    const double locusUsedDepthFraction(filteredLocusDepth * locusDepthFactor);


    if (isRNA)
    {
        double genotype(2.0);
        if (locus.is_het(sampleIndex) or locus.is_hetalt(sampleIndex)) genotype = 1.0;
        features.set(RNA_SNV_SCORING_FEATURES::GT, (genotype));

        features.set(RNA_SNV_SCORING_FEATURES::QUAL, (locus.anyVariantAlleleQuality * chromDepthFactor));
        features.set(RNA_SNV_SCORING_FEATURES::F_DP, (siteSampleInfo.n_used_calls * chromDepthFactor));
        features.set(RNA_SNV_SCORING_FEATURES::F_DPF, (siteSampleInfo.n_unused_calls * chromDepthFactor));
        features.set(RNA_SNV_SCORING_FEATURES::F_GQ, (sampleInfo.genotypeQualityPolymorphic * chromDepthFactor));
        features.set(RNA_SNV_SCORING_FEATURES::F_GQX, (sampleInfo.gqx * chromDepthFactor));

        features.set(RNA_SNV_SCORING_FEATURES::I_AvgBaseQ, (locus.avgBaseQ));
        features.set(RNA_SNV_SCORING_FEATURES::I_AvgPos, (locus.rawPos));

        features.set(RNA_SNV_SCORING_FEATURES::I_BaseQRankSum, (locus.BaseQRankSum));
        features.set(RNA_SNV_SCORING_FEATURES::I_ReadPosRankSum, (locus.ReadPosRankSum));

        features.set(RNA_SNV_SCORING_FEATURES::I_SNVHPOL, (locus.hpol));
        features.set(RNA_SNV_SCORING_FEATURES::I_SNVSB, (locus.allele.strandBias));

        features.set(RNA_SNV_SCORING_FEATURES::AD0, (r0 * chromDepthFactor));
        features.set(RNA_SNV_SCORING_FEATURES::AD1, (r1 * chromDepthFactor));

        features.set(RNA_SNV_SCORING_FEATURES::ADR, safeFrac(r0, (r0 + r1)));

        // compute any experimental features not currently used in production
        //
        if (isComputeDevelopmentFeatures)
        {
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::I_MQ, (locus.mapqRMS));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::I_MQRankSum, (locus.MQRankSum));

            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::mapqZeroFraction, mapqZeroFraction);

            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::F_DP_NORM, locusUsedDepthFraction);

            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::QUAL_NORM,
                                          (locus.anyVariantAlleleQuality * filteredLocusDepthFactor));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQX_NORM,
                                          (sampleInfo.gqx * filteredLocusDepthFactor));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQ_NORM,
                                          (sampleInfo.genotypeQualityPolymorphic * filteredLocusDepthFactor));

            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::AD0_NORM,
                                          (r0 * filteredLocusDepthFactor));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::AD1_NORM,
                                          (r1 * filteredLocusDepthFactor));

            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::QUAL_EXACT,
                                          (locus.anyVariantAlleleQuality));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQX_EXACT, (sampleInfo.gqx));
            developmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQ_EXACT, (sampleInfo.genotypeQualityPolymorphic));
        }
    }
    else
    {
        {
            double genotype(0);
            if (locus.is_hetalt(sampleIndex)) genotype = 2;
            else if (not locus.is_het(sampleIndex)) genotype = 1;
            features.set(GERMLINE_SNV_SCORING_FEATURES::GENO, genotype);
        }

        features.set(GERMLINE_SNV_SCORING_FEATURES::I_MQ, (locus.mapqRMS));
        features.set(GERMLINE_SNV_SCORING_FEATURES::I_SNVHPOL, (locus.hpol));
        features.set(GERMLINE_SNV_SCORING_FEATURES::I_SNVSB, (locus.allele.strandBias));
        features.set(GERMLINE_SNV_SCORING_FEATURES::I_MQRankSum, (locus.MQRankSum));
        features.set(GERMLINE_SNV_SCORING_FEATURES::I_ReadPosRankSum, (locus.ReadPosRankSum));

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
            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::I_BaseQRankSum, (BaseQRankSum));

            // allele bias metrics
            {
                const double allelebiaslower = cdf(boost::math::binomial(r0 + r1, 0.5), r0);
                const double allelebiasupper = cdf(boost::math::binomial(r0 + r1, 0.5), r1);

                // +1e-30 to avoid log(0) in extreme cases
                developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::ABlower, (-log(allelebiaslower + 1.e-30)));
                developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::AB,
                                              (-log(std::min(1., 2. * std::min(allelebiaslower, allelebiasupper)) + 1.e-30)));
            }

            //The average baseQ of the position of alt allele
            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::I_RawBaseQ, (locus.avgBaseQ));

            //the average position value within a read of alt allele
            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::I_RawPos, (locus.rawPos));

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
                                          (r0 * filteredLocusDepthFactor));

            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::QUAL_EXACT,
                                          (locus.anyVariantAlleleQuality));
            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQ_EXACT, (sampleInfo.genotypeQualityPolymorphic));

            developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::AD1_NORM,
                                       (r1 * filteredLocusDepthFactor));
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

    const auto& ploidy(sampleInfo.getPloidy());
    uint8_t allele0Index(0);
    uint8_t allele1Index(0);
    unsigned primaryAltAlleleIndex(0);

    if (ploidy.isDiploid())
    {
        VcfGenotypeUtil::getAlleleIndices(sampleInfo.max_gt(), allele0Index, allele1Index);
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
        VcfGenotypeUtil::getAlleleIndices(sampleInfo.max_gt(), allele0Index);
        if (allele0Index>0)
        {
            primaryAltAlleleIndex = allele0Index-1;
        }
    }

    ///TODO STREL-125 generalize to multiple alts:
    const auto& primaryAltAllele(locus.getIndelAlleles()[primaryAltAlleleIndex]);

    const unsigned confidentRefCount(
        sampleInfo.supportCounts.getCounts(true).confidentRefAlleleCount() +
        sampleInfo.supportCounts.getCounts(false).confidentRefAlleleCount());
    const unsigned confidentPrimaryAltCount(
        sampleInfo.supportCounts.getCounts(true).confidentAltAlleleCount(primaryAltAlleleIndex) +
        sampleInfo.supportCounts.getCounts(false).confidentAltAlleleCount(primaryAltAlleleIndex));

    // cdf of binomial prob of seeing no more than the number of 'allele A' reads out of A reads + B reads, given p=0.5
    // cdf of binomial prob of seeing no more than the number of 'allele B' reads out of A reads + B reads, given p=0.5
    double allelebiaslower = 0;
    double allelebiasupper = 0;
    if (ploidy.isDiploid())
    {
        const auto& fwdCounts(sampleInfo.supportCounts.getCounts(true));
        const auto& revCounts(sampleInfo.supportCounts.getCounts(false));

        if (allele0Index != allele1Index)
        {
            const double r0(fwdCounts.confidentAlleleCount(allele0Index) + revCounts.confidentAlleleCount(allele0Index));
            const double r1(fwdCounts.confidentAlleleCount(allele1Index) + revCounts.confidentAlleleCount(allele1Index));
            allelebiaslower = cdf(boost::math::binomial(r0 + r1, 0.5), r0);
            allelebiasupper = cdf(boost::math::binomial(r0 + r1, 0.5), r1);
        }
    }

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

        // +1e-30 to avoid log(0) in extreme cases
        features.set(RNA_INDEL_SCORING_FEATURES::ABlower, (-std::log(allelebiaslower + 1.e-30)));
        features.set(RNA_INDEL_SCORING_FEATURES::AB,
                     (-std::log(std::min(1., 2. * std::min(allelebiaslower, allelebiasupper)) + 1.e-30)));


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
        // code genotype as expected in the RF model:
        {
            // HET
            double genotype(0);
            if (ploidy.isDiploid())
            {
                // HOM
                if (allele0Index == allele1Index)
                {
                    genotype = 1;
                }
                else
                {
                    // HETALT
                    if ((allele0Index>0) and (allele1Index>0)) genotype = 2;
                }
            }
            features.set(GERMLINE_INDEL_SCORING_FEATURES::GENO, genotype);
        }

        features.set(GERMLINE_INDEL_SCORING_FEATURES::IDREP1, (primaryAltAllele.indelReportInfo.indel_repeat_count));
        features.set(GERMLINE_INDEL_SCORING_FEATURES::RULEN1, (primaryAltAllele.indelReportInfo.repeat_unit.length()));

        // +1e-30 to avoid log(0) in extreme cases
        features.set(GERMLINE_INDEL_SCORING_FEATURES::ABlower, (-std::log(allelebiaslower + 1.e-30)));
        features.set(GERMLINE_INDEL_SCORING_FEATURES::AB,
                     (-std::log(std::min(1., 2. * std::min(allelebiaslower, allelebiasupper)) + 1.e-30)));

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



void
GermlineDiploidIndelLocusInfo::
dump(std::ostream& os) const
{
    os << "Alleles:\n";
    for (const auto& allele : getIndelAlleles())
    {
        os << allele << "\n";
    }
}
