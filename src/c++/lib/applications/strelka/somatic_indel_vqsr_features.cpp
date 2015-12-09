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

/// \author Peter Krusche
///

#include "somatic_indel_vqsr_features.hh"
#include "somatic_call_shared.hh"
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


/**
 * Approximate indel AF from reads
 */
double
calculateIndelAF(
        const starling_indel_sample_report_info &isri
)
{
    return safeFrac(isri.n_q30_indel_reads, isri.n_q30_ref_reads + isri.n_q30_alt_reads + isri.n_q30_indel_reads);
}

/**
 * Approximate indel "other" frequency (OF) from reads
 */
double
calculateIndelOF(
        const starling_indel_sample_report_info &isri
)
{
    return safeFrac(isri.n_other_reads, isri.n_other_reads + isri.n_q30_ref_reads + isri.n_q30_alt_reads + isri.n_q30_indel_reads);
}

/**
 * Similar to
 * https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_annotator_StrandOddsRatio.php
 *
 * We adjust the counts for low coverage/low AF by adding 0.5 -- this deals with the
 * case where we don't actually observe reads on one strand.
 */
double
calculateSOR(
        const starling_indel_sample_report_info &isri
)
{

    // from
    // http://www.people.fas.harvard.edu/~mparzen/published/parzen17.pdf
    // we add .5 to each count to deal with the case of 0 / inf outcomes
    double Y1  = isri.n_q30_ref_reads_fwd + 0.5;
    double n1_minus_Y1 = isri.n_q30_indel_reads_fwd + 0.5;
    double Y2  = isri.n_q30_ref_reads_rev + 0.5;
    double n2_minus_Y2 = isri.n_q30_indel_reads_rev + 0.5;

    return log10((Y1*n1_minus_Y1)/(Y2*n2_minus_Y2));
}


/**
 * Calculate phred-scaled Fisher strand bias (p-value for the null hypothesis that
 * either REF or ALT counts are biased towards one particular strand)
 */
double
calculateFS(const starling_indel_sample_report_info & isri)
{
    return fisher_exact_test_pval_2x2(isri.n_q30_ref_reads_fwd, isri.n_q30_indel_reads_fwd,
                                      isri.n_q30_ref_reads_rev, isri.n_q30_indel_reads_rev);
}

/**
 * Calculate the p-value using a binomial test for the null hypothesis that
 * the ALT allele occurs on a particular strand only.
 */
double
calculateBSA(const starling_indel_sample_report_info & isri)
{
    return get_binomial_twosided_exact_pval(0.5, isri.n_q30_indel_reads_fwd, isri.n_q30_indel_reads) ;
}


/**
 * Calculate base-calling noise from window average set
 */
double
calculateBCNoise(const win_avg_set & was)
{
    const double filt(was.ss_filt_win.avg());
    const double used(was.ss_used_win.avg());
    const double bcnoise(safeFrac((int)filt,(int)(filt+used)));
    return bcnoise;
}

/**
 * Calculate AOR feature (log ratio between T_AF and N_AF)
 */
double
calculateAlleleFrequencyRate(const starling_indel_sample_report_info &nisri,
                             const starling_indel_sample_report_info &tisri
)
{
    const double T_AF = calculateIndelAF(tisri);
    const double N_AF = calculateIndelAF(nisri);

    return log10(std::max(T_AF, (double)0.00001) / std::max(N_AF, (double)0.0001));
}

/**
 * Calculate OD feature (log ratio between T_AF and T_OF)
 */
double
calculateTumorNoiseRate(const starling_indel_sample_report_info &tisri)
{
    const double T_AF = calculateIndelAF(tisri);
    const double T_OF = calculateIndelOF(tisri);

    return log10(std::max(T_AF, (double)0.00001) / std::max(T_OF, (double)0.0001));
}

/**
 * Calculate LAR feature (log ratio between #alt reads in tumor and #ref reads in normal)
 */
double
calculateLogAltRatio(const starling_indel_sample_report_info &nisri,
                     const starling_indel_sample_report_info &tisri)
{
    const unsigned n_ref_reads = nisri.n_q30_ref_reads;
    const unsigned t_alt_reads = tisri.n_q30_indel_reads;
    return log10(safeFrac(t_alt_reads, n_ref_reads));
}

/**
 * Calculate LOR feature (log odds ratio for  T_REF T_ALT
 *                                            N_REF N_ALT)
 */
double
calculateLogOddsRatio(const starling_indel_sample_report_info &nisri,
                      const starling_indel_sample_report_info &tisri)
{
    const double n_ref_reads = nisri.n_q30_ref_reads + .5;
    const double n_alt_reads = nisri.n_q30_indel_reads + .5;
    const double t_ref_reads = tisri.n_q30_ref_reads + .5;
    const double t_alt_reads = tisri.n_q30_indel_reads + .5;

    return log10(t_ref_reads*n_alt_reads / t_alt_reads / n_ref_reads);
}

/**
 * Calculate VQSR features and add to smod
 */
void
calculateVQSRFeatures(
    const SomaticIndelVcfInfo& siInfo,
    const win_avg_set & /* n_was */,
    const win_avg_set & /* t_was */,
    const strelka_options & /* opt */,
    const strelka_deriv_options & /* dopt */,
    strelka_shared_modifiers_indel & smod
)
{
    const somatic_indel_call::result_set& rs(siInfo.sindel.rs);
    smod.set_feature(STRELKA_INDEL_VQSR_FEATURES::QSI_NT, rs.sindel_from_ntype_qphred);
    smod.set_feature(STRELKA_INDEL_VQSR_FEATURES::ABS_T_RR, fabs(siInfo.tisri[0].readpos_ranksum.get_u_stat()));
    smod.set_feature(STRELKA_INDEL_VQSR_FEATURES::ABS_T_SOR, fabs(calculateSOR(siInfo.tisri[0])));
    smod.set_feature(STRELKA_INDEL_VQSR_FEATURES::TNR, calculateTumorNoiseRate(siInfo.tisri[0]));
    smod.set_feature(STRELKA_INDEL_VQSR_FEATURES::AFR, calculateAlleleFrequencyRate(siInfo.nisri[0], siInfo.tisri[0]));
    smod.set_feature(STRELKA_INDEL_VQSR_FEATURES::LAR, calculateLogAltRatio(siInfo.nisri[0], siInfo.tisri[0]));
    smod.set_feature(STRELKA_INDEL_VQSR_FEATURES::LOR, calculateLogOddsRatio(siInfo.nisri[0], siInfo.tisri[0]));

//    double meanChrDepth = 1;
//    auto cd = dopt.sfilter.chrom_depth.find(opt.bam_seq_name);
//    if(cd != dopt.sfilter.chrom_depth.end())
//    {
//        meanChrDepth = cd->second;
//    }

}

