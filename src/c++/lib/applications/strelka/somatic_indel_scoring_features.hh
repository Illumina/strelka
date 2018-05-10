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

#pragma once

#include "strelka_shared.hh"

#include "starling_common/AlleleReportInfo.hh"
#include "../../starling_common/LocalRegionStats.hh"
#include "somatic_call_shared.hh"
#include "somatic_indel_grid.hh"
#include "strelka_vcf_locus_info.hh"
#include "SomaticIndelVcfWriter.hh"


/// \brief Get the frequency of reads supporting the indel allele at the locus in the given sample
///
/// Note this method only considers reads which confidently support only one of the possible alleles at the locus.
/// If there are no confident reads at the locus, the method reports a frequency of 0.
///
double
getSampleIndelAlleleFrequency(
    const AlleleSampleReportInfo& isri);

/// \brief Get the frequency of reads supporting any allele other than the called indel or reference at the locus
///        in the given sample
///
/// Note this method only considers reads which confidently support only one of the possible alleles at the locus.
/// If there are no confident reads at the locus, the method reports a frequency of 0.
///
double
getSampleOtherAlleleFrequency(
    const AlleleSampleReportInfo& indelSampleReportInfo);

/// Transform input N to (N+1/N)
inline
double
makeSymmetric(const double inputRatio)
{
    assert(inputRatio > 0);
    return inputRatio + 1.0/inputRatio;
}

/// \brief Compute the allele strand odds ratio feature for the given sample.
///
/// This feature is similar to the GATK StrandOddsRatio:
/// https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_annotator_StrandOddsRatio.php
///
/// Also associated with the following publication: http://www.people.fas.harvard.edu/~mparzen/published/parzen17.pdf
///
/// Counts are adjusted for low coverage by adding a pseudocount.
///
/// \TODO Document how this feature is supposed to function and what its relationship is the various linked references.
///
double
getSampleStrandOddsRatio(
    unsigned fwdAltAlleleCount,
    unsigned revAltAlleleCount,
    unsigned fwdOtherCount,
    unsigned revOtherCount);

/// Overload to conveniently access the StrandOddsRatio directly from an AlleleSampleReportInfo object
///
/// Note this should still work for het alt cases because "ref_reads" in the AlleleSampleReportInfo typically includes
/// all alleles besides the indel in question.
///
inline
double
getSampleStrandOddsRatio(
    const AlleleSampleReportInfo& alleleSampleReportInfo)
{
    return getSampleStrandOddsRatio(
               alleleSampleReportInfo.n_confident_indel_reads_fwd,
               alleleSampleReportInfo.n_confident_indel_reads_rev,
               alleleSampleReportInfo.n_confident_ref_reads_fwd,
               alleleSampleReportInfo.n_confident_ref_reads_rev);
}

/// Calculate phred-scaled Fisher strand bias (p-value for the null hypothesis that
/// either REF or ALT counts are biased towards one particular strand)
///
double
calculateFS(const AlleleSampleReportInfo& indelSampleReportInfo);

/// Calculate the p-value using a binomial test for the null hypothesis that
/// the ALT allele occurs on a particular strand only.
///
double
calculateBSA(const AlleleSampleReportInfo& indelSampleReportInfo);

/// Calculate base-calling noise from window average set
///
double
calculateBCNoise(const LocalRegionStats& was);

/// \brief Get log ratio fo the frequency of the called indel in the tumor vs normal sample
///
double
getTumorNormalIndelAlleleLogOdds(
    const AlleleSampleReportInfo& normalIndelSampleReportInfo,
    const AlleleSampleReportInfo& tumorIndelSampleReportInfo);

/// \brief Get log ratio of the frequency of the called indel vs other non-reference alleles in the given sample
///
double
getSampleIndelNoiseLogOdds(
    const AlleleSampleReportInfo& indelSampleReportInfo);

/// Calculate LAR feature (log ratio between #alt reads in tumor and #ref reads in normal)
///
double
calculateLogAltRatio(const AlleleSampleReportInfo& nisri,
                     const AlleleSampleReportInfo& tisri);

/// \brief The log odds ratio of normalIndelAlleleCount*tumorRefAlleleCount vs normalRefAlleleCount*tumorIndelAlleleCount
///
double
getIndelAlleleCountLogOddsRatio(
    const AlleleSampleReportInfo& normalIndelSampleReportInfo,
    const AlleleSampleReportInfo& tumorIndelSampleReportInfo);


/// \brief Calculate empirical variant scoring features for candidate somatic indel
///
void
calculateScoringFeatures(
    const SomaticIndelVcfInfo& siInfo,
    const LocalRegionStats& n_was,
    const LocalRegionStats& t_was,
    const strelka_options& opt,
    strelka_shared_modifiers_indel& smod);
