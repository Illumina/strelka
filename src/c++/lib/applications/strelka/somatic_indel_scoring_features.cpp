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
/// \author Peter Krusche
///

#include "somatic_indel_scoring_features.hh"

#include "somatic_result_set.hh"
#include "somatic_indel_grid.hh"
#include "strelka_vcf_locus_info.hh"
#include "blt_util/blt_exception.hh"
#include "blt_util/io_util.hh"
#include "blt_util/qscore.hh"
#include "blt_util/fisher_exact_test.hh"
#include "blt_util/binomial_test.hh"

#include <limits>
#include <cassert>

static const double minFreq(0.0001);

static inline
double
safeFrac(const int num, const int denom)
{
    return ( (denom > 0) ? (num/static_cast<double>(denom)) : 0.);
}



double
getSampleIndelAlleleFrequency(
    const AlleleSampleReportInfo& isri)
{
    return safeFrac(isri.n_confident_indel_reads, isri.n_confident_ref_reads + isri.n_confident_alt_reads + isri.n_confident_indel_reads);
}



double
getSampleOtherAlleleFrequency(
    const AlleleSampleReportInfo& indelSampleReportInfo)
{
    return safeFrac(indelSampleReportInfo.n_other_reads, indelSampleReportInfo.n_other_reads + indelSampleReportInfo.n_confident_ref_reads + indelSampleReportInfo.n_confident_alt_reads + indelSampleReportInfo.n_confident_indel_reads);
}



double
getSampleStrandOddsRatio(
    unsigned fwdAltAlleleCount,
    unsigned revAltAlleleCount,
    unsigned fwdOtherCount,
    unsigned revOtherCount)
{
    static const double pseudocount(0.5);

    // Eq 1.1 in http://www.people.fas.harvard.edu/~mparzen/published/parzen17.pdf
    const double Y1  = fwdOtherCount + pseudocount;
    const double n1_minus_Y1 = fwdAltAlleleCount + pseudocount;
    const double Y2  = revOtherCount + pseudocount;
    const double n2_minus_Y2 = revAltAlleleCount + pseudocount;

    return (Y1*n2_minus_Y2)/(Y2*n1_minus_Y1);
}




double
calculateFS(const AlleleSampleReportInfo& indelSampleReportInfo)
{
    return fisher_exact_test_pval_2x2(indelSampleReportInfo.n_confident_ref_reads_fwd, indelSampleReportInfo.n_confident_indel_reads_fwd,
                                      indelSampleReportInfo.n_confident_ref_reads_rev, indelSampleReportInfo.n_confident_indel_reads_rev);
}



double
calculateBSA(const AlleleSampleReportInfo& indelSampleReportInfo)
{
    return get_binomial_twosided_exact_pval(0.5, indelSampleReportInfo.n_confident_indel_reads_fwd, indelSampleReportInfo.n_confident_indel_reads) ;
}



double
calculateBCNoise(const LocalRegionStats& was)
{
    const double filt(was.regionUnusedBasecallCount.avg());
    const double used(was.regionUsedBasecallCount.avg());
    const double bcnoise(safeFrac((int)filt,(int)(filt+used)));
    return bcnoise;
}



double
getTumorNormalIndelAlleleLogOdds(
    const AlleleSampleReportInfo& normalIndelSampleReportInfo,
    const AlleleSampleReportInfo& tumorIndelSampleReportInfo)
{
    const double tumorSampleIndelAlleleFrequency = getSampleIndelAlleleFrequency(tumorIndelSampleReportInfo);
    const double normalSampleIndelAlleleFrequency = getSampleIndelAlleleFrequency(normalIndelSampleReportInfo);

    return std::log(std::max(tumorSampleIndelAlleleFrequency, minFreq) / std::max(normalSampleIndelAlleleFrequency, minFreq));
}



double
getSampleIndelNoiseLogOdds(const AlleleSampleReportInfo& indelSampleReportInfo)
{
    const double indelAlleleFrequency = getSampleIndelAlleleFrequency(indelSampleReportInfo);
    const double otherAlleleFrequency = getSampleOtherAlleleFrequency(indelSampleReportInfo);

    return std::log(std::max(indelAlleleFrequency, minFreq) / std::max(otherAlleleFrequency, minFreq));
}



double
calculateLogAltRatio(
    const AlleleSampleReportInfo& normalIndelSampleReportInfo,
    const AlleleSampleReportInfo& tumorIndelSampleReportInfo)
{
    const unsigned n_ref_reads = normalIndelSampleReportInfo.n_confident_ref_reads;
    const unsigned t_alt_reads = tumorIndelSampleReportInfo.n_confident_indel_reads;
    return std::log(safeFrac(t_alt_reads, n_ref_reads));
}



double
getIndelAlleleCountLogOddsRatio(
    const AlleleSampleReportInfo& normalIndelSampleReportInfo,
    const AlleleSampleReportInfo& tumorIndelSampleReportInfo)
{
    static const double pseudoCount(0.5);
    const double normalRefCount = normalIndelSampleReportInfo.n_confident_ref_reads + pseudoCount;
    const double normalAltCount = normalIndelSampleReportInfo.n_confident_indel_reads + pseudoCount;
    const double tumorRefCount = tumorIndelSampleReportInfo.n_confident_ref_reads + pseudoCount;
    const double tumorAltCount = tumorIndelSampleReportInfo.n_confident_indel_reads + pseudoCount;

    return std::log((tumorRefCount*normalAltCount) / (tumorAltCount*normalRefCount));
}



void
calculateScoringFeatures(
    const SomaticIndelVcfInfo& siInfo,
    const LocalRegionStats& n_was,
    const LocalRegionStats& t_was,
    const strelka_options& opt,
    strelka_shared_modifiers_indel& smod)
{
    const indel_result_set& rs(siInfo.sindel.rs);

    {
        const int from_ref_qphred((rs.ntype == NTYPE::REF) ? rs.from_ntype_qphred : 0 );
        smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::SomaticIndelQualityAndHomRefGermlineGenotype, from_ref_qphred);
    }

    const double tumorSampleReadPosRankSum(siInfo.tisri[0].readpos_ranksum.get_z_stat());
    smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::TumorSampleReadPosRankSum, tumorSampleReadPosRankSum);
    smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::TumorSampleLogSymmetricStrandOddsRatio, std::log(makeSymmetric(getSampleStrandOddsRatio (siInfo.tisri[0]))));

    if (siInfo.indelReportInfo.isRepeatUnit())
    {
        smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::IndelRepeatCount, siInfo.indelReportInfo.indelRepeatCount);
        smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::RefRepeatCount, siInfo.indelReportInfo.refRepeatCount);
        smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::RepeatUnitLength, siInfo.indelReportInfo.repeatUnitLength);
    }
    else
    {
        // For complex ("SWAP") indels there will typically not be a repeat unit.
        // By making the following features default to 1 rather than 0 it becomes
        // harder for EVS to separate out and discriminate against complex indels,
        // which may be labelled incorrectly in the training data.
        smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::IndelRepeatCount, 1);
        smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::RefRepeatCount, 1);
        smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::RepeatUnitLength, 1);
    }
    smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::InterruptedHomopolymerLength, siInfo.indelReportInfo.interruptedHomopolymerLength);
    smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::TumorSampleIndelNoiseLogOdds, getSampleIndelNoiseLogOdds(siInfo.tisri[0]));
    smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::TumorNormalIndelAlleleLogOdds,
                      getTumorNormalIndelAlleleLogOdds(siInfo.nisri[0], siInfo.tisri[0]));
    smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::AlleleCountLogOddsRatio,
                      getIndelAlleleCountLogOddsRatio(siInfo.nisri[0], siInfo.tisri[0]));

    // this is how we could get the mean chromosome depth here if we needed to.
    //    double meanChrDepth = 1siInfo.iri.repeat_unit_length;
    //    auto cd = dopt.sfilter.chrom_depth.find(opt.bam_seq_name);
    //    if(cd != dopt.sfilter.chrom_depth.end())
    //    {
    //        meanChrDepth = cd->second;
    //    }

    if (opt.isReportEVSFeatures)
    {
        // features not used in the current EVS model but feature candidates/exploratory for new EVS models:
        smod.dfeatures.set(SOMATIC_INDEL_SCORING_DEVELOPMENT_FEATURES::N_AF,
                           getSampleIndelAlleleFrequency(siInfo.nisri[0]));
        smod.dfeatures.set(SOMATIC_INDEL_SCORING_DEVELOPMENT_FEATURES::T_AF,
                           getSampleIndelAlleleFrequency(siInfo.tisri[0]));

        smod.dfeatures.set(SOMATIC_INDEL_SCORING_DEVELOPMENT_FEATURES::N_OF,
                           getSampleOtherAlleleFrequency(siInfo.nisri[0]));
        smod.dfeatures.set(SOMATIC_INDEL_SCORING_DEVELOPMENT_FEATURES::T_OF,
                           getSampleOtherAlleleFrequency(siInfo.tisri[0]));

        smod.dfeatures.set(SOMATIC_INDEL_SCORING_DEVELOPMENT_FEATURES::N_BCN, calculateBCNoise(n_was));
        smod.dfeatures.set(SOMATIC_INDEL_SCORING_DEVELOPMENT_FEATURES::T_BCN, calculateBCNoise(t_was));

        smod.dfeatures.set(SOMATIC_INDEL_SCORING_DEVELOPMENT_FEATURES::TumorSampleAbsReadPosRankSum, fabs(tumorSampleReadPosRankSum));
        smod.dfeatures.set(SOMATIC_INDEL_SCORING_DEVELOPMENT_FEATURES::TumorSampleLogStrandOddsRatio, std::log(getSampleStrandOddsRatio(siInfo.tisri[0])));
        smod.dfeatures.set(SOMATIC_INDEL_SCORING_DEVELOPMENT_FEATURES::TumorSampleAbsLogStrandOddsRatio, fabs(std::log(getSampleStrandOddsRatio(siInfo.tisri[0]))));

    }
}
