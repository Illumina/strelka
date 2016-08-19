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
        os << indelSampleInfo.reportInfo << "\n";
        os << "sp size: " << indelSampleInfo.sitePloidy.size() << "\n";
    }

    return os;
}



void
GermlineDiploidSiteLocusInfo::
computeEmpiricalScoringFeatures(
    const bool isRNA,
    const bool isUniformDepthExpected,
    const bool isComputeDevelopmentFeatures,
    const double chromDepth)
{
    const double chromDepthFactor(safeFrac(1, chromDepth));

    ///TODO STREL-125 generalize to multi-sample
    const auto& firstSampleInfo(getSample(0));

    const double filteredLocusDepth(n_used_calls);
    const double locusDepth(mapqCount);

    const double filteredLocusDepthFactor(safeFrac(1, filteredLocusDepth));
    const double locusDepthFactor(safeFrac(1, locusDepth));

    // get the alt base id (choose second in case of an alt het....)
    unsigned altBase(N_BASE);
    for (unsigned b(0); b < N_BASE; ++b)
    {
        if (b == dgt.ref_gt) continue;
        if (DIGT::expect2(b, allele.max_gt))
        {
            altBase = b;
        }
    }
    assert(altBase != N_BASE);

    const unsigned r0 = alleleObservationCounts(dgt.ref_gt);
    const unsigned r1 = alleleObservationCounts(altBase);

    const double mapqZeroFraction(safeFrac(mapqZeroCount, mapqCount));

    const double locusUsedDepthFraction(filteredLocusDepth * locusDepthFactor);


    if (isRNA)
    {
        double genotype(2.0);
        if (is_het() or is_hetalt()) genotype = 1.0;
        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::GT, (genotype));

        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::QUAL, (dgt.genome.snp_qphred * chromDepthFactor));
        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::F_DP, (n_used_calls * chromDepthFactor));
        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::F_DPF, (n_unused_calls * chromDepthFactor));
        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::F_GQ, (firstSampleInfo.genotypeQualityPolymorphic * chromDepthFactor));
        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::F_GQX, (firstSampleInfo.gqx * chromDepthFactor));

        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::I_AvgBaseQ, (avgBaseQ));
        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::I_AvgPos, (rawPos));

        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::I_BaseQRankSum, (BaseQRankSum));
        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::I_ReadPosRankSum, (ReadPosRankSum));

        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::I_SNVHPOL, (hpol));
        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::I_SNVSB, (allele.strand_bias));

        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::AD0, (r0 * chromDepthFactor));
        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::AD1, (r1 * chromDepthFactor));

        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::ADR, safeFrac(r0, (r0 + r1)));

        // compute any experimental features not currently used in production
        //
        if (isComputeDevelopmentFeatures)
        {
            EVSDevelopmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::I_MQ, (mapqRMS));
            EVSDevelopmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::I_MQRankSum, (MQRankSum));

            EVSDevelopmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::mapqZeroFraction, mapqZeroFraction);

            EVSDevelopmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::F_DP_NORM, locusUsedDepthFraction);

            EVSDevelopmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::QUAL_NORM,
                                          (dgt.genome.snp_qphred * filteredLocusDepthFactor));
            EVSDevelopmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQX_NORM,
                                          (firstSampleInfo.gqx * filteredLocusDepthFactor));
            EVSDevelopmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQ_NORM,
                                          (firstSampleInfo.genotypeQualityPolymorphic * filteredLocusDepthFactor));

            EVSDevelopmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::AD0_NORM,
                                          (r0 * filteredLocusDepthFactor));
            EVSDevelopmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::AD1_NORM,
                                          (r1 * filteredLocusDepthFactor));

            EVSDevelopmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::QUAL_EXACT,
                                          (dgt.genome.snp_qphred));
            EVSDevelopmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQX_EXACT, (firstSampleInfo.gqx));
            EVSDevelopmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQ_EXACT, (firstSampleInfo.genotypeQualityPolymorphic));
        }
    }
    else
    {
        {
            double genotype(0);
            if (is_hetalt()) genotype = 2;
            else if (not is_het()) genotype = 1;
            EVSFeatures.set(GERMLINE_SNV_SCORING_FEATURES::GENO, genotype);
        }

        EVSFeatures.set(GERMLINE_SNV_SCORING_FEATURES::I_MQ, (mapqRMS));
        EVSFeatures.set(GERMLINE_SNV_SCORING_FEATURES::I_SNVHPOL, (hpol));
        EVSFeatures.set(GERMLINE_SNV_SCORING_FEATURES::I_SNVSB, (allele.strand_bias));
        EVSFeatures.set(GERMLINE_SNV_SCORING_FEATURES::I_MQRankSum, (MQRankSum));
        EVSFeatures.set(GERMLINE_SNV_SCORING_FEATURES::I_ReadPosRankSum, (ReadPosRankSum));

        // how surprising is the depth relative to expect? This is the only value will be modified for exome/targeted runs
        /// TODO: convert this to pvalue based on Poisson distro?
        double relativeLocusDepth(1.);
        if (isUniformDepthExpected)
        {
            relativeLocusDepth = (locusDepth * chromDepthFactor);
        }

        EVSFeatures.set(GERMLINE_SNV_SCORING_FEATURES::TDP_NORM, relativeLocusDepth);

        // how noisy is the locus?
        EVSFeatures.set(GERMLINE_SNV_SCORING_FEATURES::F_DP_NORM, locusUsedDepthFraction);

        EVSFeatures.set(GERMLINE_SNV_SCORING_FEATURES::F_GQX_EXACT, (firstSampleInfo.gqx));

        // compute any experimental features not currently used in production
        //
        if (isComputeDevelopmentFeatures)
        {

            // BaseQRankSum
            EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::I_BaseQRankSum, (BaseQRankSum));

            // allele bias metrics
            {
                const double allelebiaslower = cdf(boost::math::binomial(r0 + r1, 0.5), r0);
                const double allelebiasupper = cdf(boost::math::binomial(r0 + r1, 0.5), r1);

                // +1e-30 to avoid log(0) in extreme cases
                EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::ABlower, (-log(allelebiaslower + 1.e-30)));
                EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::AB,
                                              (-log(std::min(1., 2. * std::min(allelebiaslower, allelebiasupper)) + 1.e-30)));
            }

            //The average baseQ of the position of alt allele
            EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::I_RawBaseQ, (avgBaseQ));

            //the average position value within a read of alt allele
            EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::I_RawPos, (rawPos));

            // hom unrelable are the read mappings near this locus?
            EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::mapqZeroFraction, mapqZeroFraction);

            // renormalized features intended to replace the corresponding production feature
            //
            EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::QUAL_NORM,
                                          (dgt.genome.snp_qphred * filteredLocusDepthFactor));
            EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQX_NORM,
                                          (firstSampleInfo.gqx * filteredLocusDepthFactor));
            EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQ_NORM,
                                          (firstSampleInfo.genotypeQualityPolymorphic * filteredLocusDepthFactor));

            EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::AD0_NORM,
                                          (r0 * filteredLocusDepthFactor));

            EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::QUAL_EXACT,
                                          (dgt.genome.snp_qphred));
       
            EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQ_EXACT, (firstSampleInfo.genotypeQualityPolymorphic));

            EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::AD1_NORM,
                                          (r1 * filteredLocusDepthFactor));
        }
    }
}



std::ostream&
operator<<(std::ostream& os,
           const GermlineDiploidSiteLocusInfo& si)
{
    os << "pos: " << (si.pos+1) << " " << si.get_gt();
    return os;
}



void
GermlineDiploidIndelLocusInfo::
computeEmpiricalScoringFeatures(
    const bool isRNA,
    const bool isUniformDepthExpected,
    const bool isComputeDevelopmentFeatures,
    const double chromDepth)
{
    ///TODO STREL-125 generalize to multi-sample
    const auto& firstSampleInfo(getSample(0));
    const auto& firstIndelSampleInfo(getIndelSample(0));

    ///TODO STREL-125 generalize to multiple alts:
    const auto& firstAltAllele(getIndelAlleles().front());

    const auto& sampleReportInfo(firstIndelSampleInfo.reportInfo);

    const double filteredLocusDepth(sampleReportInfo.tier1Depth);
    const double locusDepth(sampleReportInfo.mapqTracker.count);
    const double confidentDepth(sampleReportInfo.total_confident_reads());

    const double chromDepthFactor(safeFrac(1,chromDepth));
    const double filteredLocusDepthFactor(safeFrac(1,filteredLocusDepth));
    const double locusDepthFactor(safeFrac(1,locusDepth));
    const double confidentDepthFactor(safeFrac(1,confidentDepth));

    const bool isDiploid(firstSampleInfo.getPloidy().isDiploid());
    uint8_t allele0Index(0);
    uint8_t allele1Index(0);

    if (isDiploid)
    {
        VcfGenotypeUtil::getAlleleIndices(firstSampleInfo.max_gt(), allele0Index, allele1Index);
    }

    // cdf of binomial prob of seeing no more than the number of 'allele A' reads out of A reads + B reads, given p=0.5
    // cdf of binomial prob of seeing no more than the number of 'allele B' reads out of A reads + B reads, given p=0.5
    double allelebiaslower = 0;
    double allelebiasupper = 0;
    if (isDiploid)
    {
        const auto& fwdCounts(firstSampleInfo.supportCounts.getCounts(true));
        const auto& revCounts(firstSampleInfo.supportCounts.getCounts(false));

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
        features.set(RNA_INDEL_SCORING_FEATURES::QUAL, (anyVariantAlleleQuality * chromDepthFactor));
        features.set(RNA_INDEL_SCORING_FEATURES::F_GQX, (firstSampleInfo.gqx * chromDepthFactor));
        features.set(RNA_INDEL_SCORING_FEATURES::REFREP1, (firstAltAllele.indelReportInfo.ref_repeat_count));
        features.set(RNA_INDEL_SCORING_FEATURES::IDREP1, (firstAltAllele.indelReportInfo.indel_repeat_count));
        features.set(RNA_INDEL_SCORING_FEATURES::RULEN1, (firstAltAllele.indelReportInfo.repeat_unit.length()));
        features.set(RNA_INDEL_SCORING_FEATURES::AD0,
                     (sampleReportInfo.n_confident_ref_reads * chromDepthFactor));
        features.set(RNA_INDEL_SCORING_FEATURES::AD1,
                     (sampleReportInfo.n_confident_indel_reads * chromDepthFactor));
        features.set(RNA_INDEL_SCORING_FEATURES::AD2,
                     (sampleReportInfo.n_confident_alt_reads * chromDepthFactor));
        features.set(RNA_INDEL_SCORING_FEATURES::F_DPI, (sampleReportInfo.tier1Depth * chromDepthFactor));

        // +1e-30 to avoid log(0) in extreme cases
        features.set(RNA_INDEL_SCORING_FEATURES::ABlower, (-std::log(allelebiaslower + 1.e-30)));
        features.set(RNA_INDEL_SCORING_FEATURES::AB,
                     (-std::log(std::min(1., 2. * std::min(allelebiaslower, allelebiasupper)) + 1.e-30)));


        // compute any experimental features not currently used in production
        if (isComputeDevelopmentFeatures)
        {
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQ, (firstSampleInfo.genotypeQualityPolymorphic * chromDepthFactor));

            // how unreliable are the read mappings near this locus?
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_MQ,
                                    (sampleReportInfo.mapqTracker.getRMS()));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::mapqZeroFraction,
                                    (sampleReportInfo.mapqTracker.getZeroFrac()));

            // how noisy is the locus?
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_DPI_NORM,
                                    (filteredLocusDepth * locusDepthFactor));


            // all of the features below are simply renormalized replacements of the current production feature set
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::QUAL_NORM,
                                    (anyVariantAlleleQuality * filteredLocusDepthFactor));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQX_NORM,
                                    (firstSampleInfo.gqx * filteredLocusDepthFactor));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQ_NORM,
                                    (firstSampleInfo.genotypeQualityPolymorphic * filteredLocusDepthFactor));

            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::AD0_NORM,
                                    (sampleReportInfo.n_confident_ref_reads * confidentDepthFactor));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::AD1_NORM,
                                    (sampleReportInfo.n_confident_indel_reads * confidentDepthFactor));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::AD2_NORM,
                                    (sampleReportInfo.n_confident_alt_reads * confidentDepthFactor));

            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::QUAL_EXACT, (anyVariantAlleleQuality));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQX_EXACT, (firstSampleInfo.gqx));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQ_EXACT, (firstSampleInfo.genotypeQualityPolymorphic));
        }
    }
    else
    {
        // code genotype as expected in the RF model:
        {
            // HET
            double genotype(0);
            if (isDiploid)
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

        features.set(GERMLINE_INDEL_SCORING_FEATURES::IDREP1, (firstAltAllele.indelReportInfo.indel_repeat_count));
        features.set(GERMLINE_INDEL_SCORING_FEATURES::RULEN1, (firstAltAllele.indelReportInfo.repeat_unit.length()));

        // +1e-30 to avoid log(0) in extreme cases
        features.set(GERMLINE_INDEL_SCORING_FEATURES::ABlower, (-std::log(allelebiaslower + 1.e-30)));
        features.set(GERMLINE_INDEL_SCORING_FEATURES::AB,
                     (-std::log(std::min(1., 2. * std::min(allelebiaslower, allelebiasupper)) + 1.e-30)));

        // how unreliable are the read mappings near this locus?
        features.set(GERMLINE_INDEL_SCORING_FEATURES::F_MQ,
                     (sampleReportInfo.mapqTracker.getRMS()));

        // how surprising is the depth relative to expect? This is the only value will be modified for exome/targeted runs
        //
        /// TODO: convert this to pvalue based on Poisson distro?
        double relativeLocusDepth(1.);
        if (isUniformDepthExpected)
        {
            relativeLocusDepth = (locusDepth * chromDepthFactor);
        }
        features.set(GERMLINE_INDEL_SCORING_FEATURES::TDP_NORM, relativeLocusDepth);

        features.set(GERMLINE_INDEL_SCORING_FEATURES::AD1_NORM,
                     (sampleReportInfo.n_confident_indel_reads * confidentDepthFactor));

        features.set(GERMLINE_INDEL_SCORING_FEATURES::F_GQX_EXACT, (firstSampleInfo.gqx));

        // compute any experimental features not currently used in production
        if (isComputeDevelopmentFeatures)
        {
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::REFREP1, (firstAltAllele.indelReportInfo.ref_repeat_count));

            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::mapqZeroFraction,
                                    (sampleReportInfo.mapqTracker.getZeroFrac()));

            // how noisy is the locus?
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::F_DPI_NORM,
                                    (filteredLocusDepth * locusDepthFactor));

            // all of the features below are simply renormalized replacements of the current production feature set
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::QUAL_NORM,
                                    (anyVariantAlleleQuality * filteredLocusDepthFactor));
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQX_NORM,
                                    (firstSampleInfo.gqx * filteredLocusDepthFactor));
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQ_NORM,
                                    (firstSampleInfo.genotypeQualityPolymorphic * filteredLocusDepthFactor));

            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::AD0_NORM,
                                    (sampleReportInfo.n_confident_ref_reads * confidentDepthFactor));
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::AD2_NORM,
                                    (sampleReportInfo.n_confident_alt_reads * confidentDepthFactor));

            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::QUAL_EXACT, (anyVariantAlleleQuality));
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQ_EXACT, (firstSampleInfo.genotypeQualityPolymorphic));
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
