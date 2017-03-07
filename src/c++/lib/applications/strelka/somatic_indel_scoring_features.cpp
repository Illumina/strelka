// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
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


static inline
double
safeFrac(const int num, const int denom)
{
    return ( (denom > 0) ? (num/static_cast<double>(denom)) : 0.);
}



double
calculateIndelAF(
    const AlleleSampleReportInfo& isri)
{
    return safeFrac(isri.n_confident_indel_reads, isri.n_confident_ref_reads + isri.n_confident_alt_reads + isri.n_confident_indel_reads);
}



double
calculateIndelOF(
    const AlleleSampleReportInfo& isri)
{
    return safeFrac(isri.n_other_reads, isri.n_other_reads + isri.n_confident_ref_reads + isri.n_confident_alt_reads + isri.n_confident_indel_reads);
}



double
calculateSOR(
    const AlleleSampleReportInfo& isri)
{
    // from
    // http://www.people.fas.harvard.edu/~mparzen/published/parzen17.pdf
    // we add .5 to each count to deal with the case of 0 / inf outcomes
    double Y1  = isri.n_confident_ref_reads_fwd + 0.5;
    double n1_minus_Y1 = isri.n_confident_indel_reads_fwd + 0.5;
    double Y2  = isri.n_confident_ref_reads_rev + 0.5;
    double n2_minus_Y2 = isri.n_confident_indel_reads_rev + 0.5;

    return log10((Y1*n1_minus_Y1)/(Y2*n2_minus_Y2));
}



double
calculateFS(const AlleleSampleReportInfo& isri)
{
    return fisher_exact_test_pval_2x2(isri.n_confident_ref_reads_fwd, isri.n_confident_indel_reads_fwd,
                                      isri.n_confident_ref_reads_rev, isri.n_confident_indel_reads_rev);
}



double
calculateBSA(const AlleleSampleReportInfo& isri)
{
    return get_binomial_twosided_exact_pval(0.5, isri.n_confident_indel_reads_fwd, isri.n_confident_indel_reads) ;
}



double
calculateBCNoise(const win_avg_set& was)
{
    const double filt(was.ss_filt_win.avg());
    const double used(was.ss_used_win.avg());
    const double bcnoise(safeFrac((int)filt,(int)(filt+used)));
    return bcnoise;
}



double
calculateAlleleFrequencyRate(
    const AlleleSampleReportInfo& normalIndelSampleReportInfo,
    const AlleleSampleReportInfo& tumorIndelSampleReportInfo)
{
    const double T_AF = calculateIndelAF(tumorIndelSampleReportInfo);
    const double N_AF = calculateIndelAF(normalIndelSampleReportInfo);

    return log10(std::max(T_AF, (double)0.00001) / std::max(N_AF, (double)0.0001));
}



double
calculateTumorNoiseRate(const AlleleSampleReportInfo& tumorIndelSampleReportInfo)
{
    const double T_AF = calculateIndelAF(tumorIndelSampleReportInfo);
    const double T_OF = calculateIndelOF(tumorIndelSampleReportInfo);

    return log10(std::max(T_AF, (double)0.00001) / std::max(T_OF, (double)0.0001));
}



double
calculateLogAltRatio(
    const AlleleSampleReportInfo& normalIndelSampleReportInfo,
    const AlleleSampleReportInfo& tumorIndelSampleReportInfo)
{
    const unsigned n_ref_reads = normalIndelSampleReportInfo.n_confident_ref_reads;
    const unsigned t_alt_reads = tumorIndelSampleReportInfo.n_confident_indel_reads;
    return log10(safeFrac(t_alt_reads, n_ref_reads));
}



double
calculateLogOddsRatio(
    const AlleleSampleReportInfo& normalIndelSampleReportInfo,
    const AlleleSampleReportInfo& tumorIndelSampleReportInfo)
{
    const double n_ref_reads = normalIndelSampleReportInfo.n_confident_ref_reads + .5;
    const double n_alt_reads = normalIndelSampleReportInfo.n_confident_indel_reads + .5;
    const double t_ref_reads = tumorIndelSampleReportInfo.n_confident_ref_reads + .5;
    const double t_alt_reads = tumorIndelSampleReportInfo.n_confident_indel_reads + .5;

    return log10(t_ref_reads*n_alt_reads / t_alt_reads / n_ref_reads);
}



void
calculateScoringFeatures(
    const SomaticIndelVcfInfo& siInfo,
    const win_avg_set& n_was,
    const win_avg_set& t_was,
    const strelka_options& opt,
    const strelka_deriv_options& /* dopt */,
    strelka_shared_modifiers_indel& smod)
{
    const indel_result_set& rs(siInfo.sindel.rs);
    smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::QSI_NT, rs.from_ntype_qphred);
    smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::ABS_T_RR, fabs(siInfo.tisri[0].readpos_ranksum.get_z_stat()));
    smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::ABS_T_SOR, fabs(calculateSOR(siInfo.tisri[0])));

    if (siInfo.indelReportInfo.is_repeat_unit())
    {
        smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::IC, siInfo.indelReportInfo.indel_repeat_count);
        smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::RC, siInfo.indelReportInfo.ref_repeat_count);
        smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::RU_LEN, siInfo.indelReportInfo.repeat_unit_length);
    }
    else
    {
        smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::IC, 0);
        smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::RC, 0);
        smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::RU_LEN, 0);
    }
    smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::IHP, siInfo.indelReportInfo.ihpol);
    smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::TNR, calculateTumorNoiseRate(siInfo.tisri[0]));
    smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::AFR, calculateAlleleFrequencyRate(siInfo.nisri[0], siInfo.tisri[0]));
    smod.features.set(SOMATIC_INDEL_SCORING_FEATURES::LOR, calculateLogOddsRatio(siInfo.nisri[0], siInfo.tisri[0]));

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
        smod.dfeatures.set(SOMATIC_INDEL_SCORING_DEVELOPMENT_FEATURES::N_AF, calculateIndelAF(siInfo.nisri[0]));
        smod.dfeatures.set(SOMATIC_INDEL_SCORING_DEVELOPMENT_FEATURES::T_AF, calculateIndelAF(siInfo.tisri[0]));

        smod.dfeatures.set(SOMATIC_INDEL_SCORING_DEVELOPMENT_FEATURES::N_OF, calculateIndelOF(siInfo.nisri[0]));
        smod.dfeatures.set(SOMATIC_INDEL_SCORING_DEVELOPMENT_FEATURES::T_OF, calculateIndelOF(siInfo.tisri[0]));

        smod.dfeatures.set(SOMATIC_INDEL_SCORING_DEVELOPMENT_FEATURES::N_BCN, calculateBCNoise(n_was));
        smod.dfeatures.set(SOMATIC_INDEL_SCORING_DEVELOPMENT_FEATURES::T_BCN, calculateBCNoise(t_was));
    }
}

