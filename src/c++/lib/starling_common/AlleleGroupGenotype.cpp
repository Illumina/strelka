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

#include "AlleleGroupGenotype.hh"

#include "starling_indel_call_pprob_digt.hh"
#include "blt_util/math_util.hh"
#include "blt_util/prob_util.hh"
#include "blt_util/qscore.hh"



static
double
integrate_out_sites(
    const starling_base_deriv_options& dopt,
    const uint16_t nsite,
    const double p_on_site,
    const bool is_tier2_pass)
{
    return log_sum((p_on_site + dopt.site_lnprior),
                   (dopt.get_nonsite_path_lnp(is_tier2_pass,nsite) + dopt.nonsite_lnprior));
}



static
void
accumulateLogLhoodFromReadObservation(
    const starling_base_deriv_options& dopt,
    const starling_sample_options& sampleOptions,
    const uint16_t nsite,
    const uint16_t read_length,
    const bool is3AlleleModel,
    const double refAllele_lnp,
    const double variantAllele0_lnp,
    const double variantAllele1_lnp,
    const IndelKey& variantAllele0Key,
    const IndelKey& variantAllele1Key,
    double* logLhood)
{
    static const bool isTier2Pass(false);

    logLhood[AG_GENOTYPE::HOMREF] += integrate_out_sites(dopt,nsite,refAllele_lnp, isTier2Pass);
    logLhood[AG_GENOTYPE::HOM0] += integrate_out_sites(dopt,nsite,variantAllele0_lnp, isTier2Pass);

    static const double loghalf(std::log(0.5));
    double logHet0RefPrior(loghalf);
    double logHet0Allele0Prior(loghalf);
    {
        static const double hetAlleleRatio(0.5);
        get_het_observed_allele_ratio(read_length, sampleOptions.min_read_bp_flank,
                                      variantAllele0Key, hetAlleleRatio, logHet0RefPrior, logHet0Allele0Prior);
    }

    const double het0_lnp(log_sum(refAllele_lnp+logHet0RefPrior,variantAllele0_lnp+logHet0Allele0Prior));
    logLhood[AG_GENOTYPE::HET0] += integrate_out_sites(dopt,nsite,het0_lnp, isTier2Pass);

    if (is3AlleleModel)
    {
        logLhood[AG_GENOTYPE::HOM1] += integrate_out_sites(dopt,nsite,variantAllele1_lnp, isTier2Pass);

        double logHet1RefPrior(loghalf);
        double logHet1Allele1Prior(loghalf);
        {
            static const double hetAlleleRatio(0.5);
            get_het_observed_allele_ratio(read_length, sampleOptions.min_read_bp_flank,
                                          variantAllele1Key, hetAlleleRatio, logHet1RefPrior, logHet1Allele1Prior);
        }

        const double het1_lnp(log_sum(refAllele_lnp+logHet1RefPrior,variantAllele1_lnp+logHet1Allele1Prior));
        logLhood[AG_GENOTYPE::HET1] += integrate_out_sites(dopt,nsite,het1_lnp, isTier2Pass);

        /// approximate the expected allele ratio in this case
        const double normalizeHetRatio(log_sum(logHet0Allele0Prior, logHet1Allele1Prior));
        const double het01_lnp(log_sum(variantAllele0_lnp+(logHet0Allele0Prior-normalizeHetRatio),variantAllele1_lnp+(logHet1Allele1Prior-normalizeHetRatio)));
        logLhood[AG_GENOTYPE::HET01] += integrate_out_sites(dopt,nsite,het01_lnp, isTier2Pass);
    }
}



static
void
logLhoodToLocusGenotype(
    const ContextGenotypePriors& genotypePriors,
    const unsigned gtCount,
    const bool is3AlleleModel,
    const bool isHaploid,
    double* logLhood,
    AlleleGroupGenotype& locusGenotype)
{
    // mult by prior distro to get unnormalized pprob:
    const double* genotypeLogPrior(genotypePriors.get2Allele(isHaploid));
    if (is3AlleleModel)
    {
        genotypeLogPrior = (genotypePriors.get3Allele(isHaploid));
    }

    for (unsigned gt(0); gt<gtCount; ++gt)
    {
        locusGenotype.posteriorProb[gt] = logLhood[gt] + genotypeLogPrior[gt];
    }

    normalize_ln_distro(std::begin(locusGenotype.posteriorProb),std::begin(locusGenotype.posteriorProb)+gtCount,locusGenotype.maxGenotypeIndex);

    locusGenotype.anyVariantAlleleQuality=error_prob_to_qphred(locusGenotype.posteriorProb[AG_GENOTYPE::HOMREF]);
    locusGenotype.genotypeQuality=error_prob_to_qphred(prob_comp(locusGenotype.posteriorProb,locusGenotype.posteriorProb+gtCount,locusGenotype.maxGenotypeIndex));


    // set phredLoghood:
    {
        unsigned maxIndex(0);
        for (unsigned gt(1); gt<gtCount; ++gt)
        {
            if (logLhood[gt] > logLhood[maxIndex]) maxIndex = gt;
        }
        for (unsigned gt(0); gt<gtCount; ++gt)
        {
            // don't enforce maxQ at this point, b/c we're going to possibly select down from this list:
            locusGenotype.phredLoghood[gt] = ln_error_prob_to_qphred(logLhood[gt]-logLhood[maxIndex]);
        }
    }

    // add new poly calls (this trashes lhood):
    {
        // mult by prior distro to get unnormalized pprob:
        const double* genotypeLogPriorPolymorphic(genotypePriors.get2AllelePolymorphic(isHaploid));
        if (is3AlleleModel)
        {
            genotypeLogPriorPolymorphic = (genotypePriors.get3AllelePolymorphic(isHaploid));
        }

        for (unsigned gt(0); gt<gtCount; ++gt)
        {
            logLhood[gt] += genotypeLogPriorPolymorphic[gt];
        }
    }
    normalize_ln_distro(logLhood,logLhood+gtCount,locusGenotype.maxGenotypeIndexPolymorphic);
    locusGenotype.genotypeQualityPolymorphic=error_prob_to_qphred(prob_comp(logLhood,logLhood+gtCount,locusGenotype.maxGenotypeIndexPolymorphic));
}



///
/// \param alleleLhood[inout] allele likelihoods, note that these are altered to include nsite and
///                           be normalized in this func, 0 index is reference
/// \param locusReadStats
///
void
updateSupportingReadStats(
    const starling_base_deriv_options& dopt,
    const double readSupportTheshold,
    const uint16_t nsite,
    const bool isFwdStrand,
    std::vector<double>& alleleLoglhoods,
    LocusSupportingReadStats& locusReadStats)
{
    static const bool isTier2Pass(false);
    for (auto& alleleHood : alleleLoglhoods)
    {
        alleleHood = integrate_out_sites(dopt, nsite, alleleHood, isTier2Pass);
    }
    unsigned maxIndex(0);
    normalize_ln_distro(alleleLoglhoods.begin(),alleleLoglhoods.end(),maxIndex);

    bool isConfidentAlleleFound(false);
    const unsigned fullAlleleCount(alleleLoglhoods.size());
    for (unsigned alleleIndex(0); alleleIndex<fullAlleleCount; ++alleleIndex)
    {
        if (alleleLoglhoods[alleleIndex] < readSupportTheshold) continue;
        locusReadStats.getCounts(isFwdStrand).incrementAlleleCount(alleleIndex);
        isConfidentAlleleFound=true;
        break;
    }

    if (not isConfidentAlleleFound)
    {
        locusReadStats.getCounts(isFwdStrand).nonConfidentCount++;
    }
}



void
getVariantAlleleGroupGenotypeLhoods(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const starling_sample_options& sampleOptions,
    const reference_contig_segment& ref,
    const unsigned groupLocusPloidy,
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    AlleleGroupGenotype& locusGenotype,
    LocusSupportingReadStats& locusReadStats)
{
    assert(groupLocusPloidy>0u);
    assert(alleleGroup.size()<=groupLocusPloidy);
    assert(groupLocusPloidy<3u);

    const uint8_t nonRefAlleleCount(alleleGroup.size());
    const uint8_t fullAlleleCount(nonRefAlleleCount+1);

    GenotypeInfo ginfo(groupLocusPloidy, fullAlleleCount);
    const unsigned gtCount(ginfo.genotypeCount());
    assert(gtCount<=AlleleGroupGenotype::MAX_GENOTYPE_COUNT);
    double logLhood[AlleleGroupGenotype::MAX_GENOTYPE_COUNT];
    std::fill(logLhood,logLhood+gtCount,0.);
    const bool is3AlleleModel(fullAlleleCount == 3);

    //---------------------------------------------------
    // get patternRepeatCount, punt on figuring out how to do this for more than one alt allele
    //
    /// TODO set PatternRepeatCount for multiple alts
    //
    unsigned patternRepeatCount(1);
    if (nonRefAlleleCount==1)
    {
        const unsigned nonRefAllele0Index(0);
        const IndelKey& allele0Key(alleleGroup.key(nonRefAllele0Index));
        const IndelData& allele0Data(alleleGroup.data(nonRefAllele0Index));

        starling_indel_report_info indelReportInfo;
        get_starling_indel_report_info(allele0Key, allele0Data, ref, indelReportInfo);
        patternRepeatCount=std::max(1u,indelReportInfo.ref_repeat_count);
    }

    const ContextGenotypePriors& genotypePriors(dopt.getIndelGenotypePriors().getContextSpecificPriorSet(patternRepeatCount));

    //---------------------------------------------------
    // get genotype lhoods
    //
    if (nonRefAlleleCount>0)
    {
        const unsigned nonRefAllele0Index(0);
        const IndelKey& allele0Key(alleleGroup.key(nonRefAllele0Index));
        const IndelData& allele0Data(alleleGroup.data(nonRefAllele0Index));
        const IndelSampleData& allele0SampleData(allele0Data.getSampleData(sampleIndex));

        const unsigned nonRefAllele1Index(is3AlleleModel ? 1 : 0);
        const IndelKey& allele1Key(alleleGroup.key(nonRefAllele1Index));
        const IndelData& allele1Data(alleleGroup.data(nonRefAllele1Index));
        const IndelSampleData& allele1SampleData(allele1Data.getSampleData(sampleIndex));

        // threshold used to generate supporting count summary (not used for GT likelihoods):
        const double readSupportTheshold(opt.readConfidentSupportThreshold.numval());
        // used for support counts only:
        std::vector<double> alleleLhood(fullAlleleCount);

        locusReadStats.setAltCount(nonRefAlleleCount);

        for (const auto& score : allele0SampleData.read_path_lnp)
        {
            const auto readIndex(score.first);
            const ReadPathScores& allele0ReadScores(score.second);

            if (!allele0ReadScores.is_tier1_read) continue;

            const ReadPathScores* allele1ReadScoresPtr(nullptr);
            if (is3AlleleModel)
            {
                /// skip this read if it's not aligned to both alleles:
                const auto allele1ReadScoresIter(allele1SampleData.read_path_lnp.find(readIndex));
                if (allele1ReadScoresIter == allele1SampleData.read_path_lnp.end()) continue;
                allele1ReadScoresPtr = (&(allele1ReadScoresIter->second));
            }

            double refAllele_lnp(allele0ReadScores.ref);
            const double variantAllele0_lnp(allele0ReadScores.indel);
            double variantAllele1_lnp(variantAllele0_lnp);

            if (is3AlleleModel)
            {
                const ReadPathScores& allele1ReadScores(*allele1ReadScoresPtr);
                refAllele_lnp = std::max(refAllele_lnp, (double) allele1ReadScores.ref);
                variantAllele1_lnp = allele1ReadScores.indel;
            }

            accumulateLogLhoodFromReadObservation(
                dopt, sampleOptions,
                allele0ReadScores.nsite, allele0ReadScores.read_length,
                is3AlleleModel,
                refAllele_lnp, variantAllele0_lnp, variantAllele1_lnp,
                allele0Key, allele1Key, logLhood);

            //---------------------------------------------------
            // get support counts for each allele
            //

            alleleLhood[0] = refAllele_lnp;
            alleleLhood[1] = variantAllele0_lnp;
            if (is3AlleleModel)
            {
                alleleLhood[2] = variantAllele1_lnp;
            }
            updateSupportingReadStats(dopt, readSupportTheshold, allele0ReadScores.nsite, allele0ReadScores.is_fwd_strand, alleleLhood, locusReadStats);
        }
    }

    const bool isHaploid(groupLocusPloidy==1);
    logLhoodToLocusGenotype(genotypePriors, gtCount, is3AlleleModel, isHaploid, logLhood, locusGenotype);
}



/// for a specific read, get the likelihood of the read conditioned on a whole group of alleles
///
/// this type of likelihood can be useful for something like a VCF <*> abstract allele type
///
/// \return the variantAlleleIndex of the most likely allele for this read
///
static
unsigned
updateFromAlleleGroup(
    const unsigned sampleId,
    const IndelKey& excludeAlleleKey,
    const unsigned targetReadIndex,
    const OrthogonalVariantAlleleCandidateGroup& variantAlleleGroup,
    double& refAllele_lnp,
    double& variantAlleleGroup_lnp)
{
    bool isAnyAlleleEligible(false);
    bool isAnyAlleleScored(false);
    unsigned maxVariantAlleleIndex(0);

    const double variantAlleleCount(variantAlleleGroup.size());
    for (unsigned variantAlleleIndex(0); variantAlleleIndex<variantAlleleCount; ++variantAlleleIndex)
    {
        const IndelKey& variantAlleleKey(variantAlleleGroup.key(variantAlleleIndex));
        if (variantAlleleKey == excludeAlleleKey) continue;

        isAnyAlleleEligible = true;

        const IndelData& id(variantAlleleGroup.data(variantAlleleIndex));
        const IndelSampleData& isd(id.getSampleData(sampleId));

        const auto AlleleReadScoresIter(isd.read_path_lnp.find(targetReadIndex));
        if (AlleleReadScoresIter == isd.read_path_lnp.end()) continue;

        isAnyAlleleScored=true;

        const ReadPathScores& AlleleReadScores(AlleleReadScoresIter->second);
        refAllele_lnp = std::max(refAllele_lnp, (double) AlleleReadScores.ref);
        if (AlleleReadScores.indel > variantAlleleGroup_lnp)
        {
            variantAlleleGroup_lnp = AlleleReadScores.indel;
            maxVariantAlleleIndex = variantAlleleIndex;
        }
    }

    if (isAnyAlleleEligible and (not isAnyAlleleScored))
    {
        variantAlleleGroup_lnp = std::max(variantAlleleGroup_lnp, refAllele_lnp);
    }

    return maxVariantAlleleIndex;
}



void
getGenotypeLhoodsForForcedOutputAllele(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const starling_sample_options& sampleOptions,
    const reference_contig_segment& ref,
    const unsigned groupLocusPloidy,
    const unsigned sampleId,
    const OrthogonalVariantAlleleCandidateGroup& variantAlleleGroup,
    const OrthogonalVariantAlleleCandidateGroup& forcedOutputAlleleGroup,
    const unsigned forcedOutputAlleleIndex,
    AlleleGroupGenotype& locusGenotype,
    LocusSupportingReadStats& locusReadStats)
{
    static const double log0(-std::numeric_limits<double>::infinity());

    assert(groupLocusPloidy>0u);
    assert(groupLocusPloidy<3u);

    // this is a fixed allele set: {ref, forced output allele, everything else}
    const uint8_t fullAlleleCount(3);
    const bool is3AlleleModel(fullAlleleCount==3);

    // TEMPORARY -- we compress out the <*> allele in final output:
    const uint8_t tmpFullAlleleCount(2);
    const uint8_t tmpNonRefAlleleCount(1);

    // threshold used to generate supporting count summary (not used for GT likelihoods):
    const double readSupportTheshold(opt.readConfidentSupportThreshold.numval());
    // used for support counts only:
    std::vector<double> alleleLhood(tmpFullAlleleCount);
    locusReadStats.setAltCount(tmpNonRefAlleleCount);

    GenotypeInfo ginfo(groupLocusPloidy,fullAlleleCount);
    const unsigned gtCount(ginfo.genotypeCount());
    assert(gtCount<=AlleleGroupGenotype::MAX_GENOTYPE_COUNT);
    double logLhood[AlleleGroupGenotype::MAX_GENOTYPE_COUNT];
    std::fill(logLhood,logLhood+gtCount,0.);

    const IndelKey& forcedOutputIndelKey(forcedOutputAlleleGroup.key(forcedOutputAlleleIndex));
    const IndelData& forcedOutputIndelData(forcedOutputAlleleGroup.data(forcedOutputAlleleIndex));
    const IndelSampleData& forcedOutputIndelSampleData(forcedOutputIndelData.getSampleData(sampleId));

    unsigned patternRepeatCount=1;
    {
        starling_indel_report_info indelReportInfo;
        get_starling_indel_report_info(forcedOutputIndelKey, forcedOutputIndelData, ref, indelReportInfo);
        patternRepeatCount=std::max(1u,indelReportInfo.ref_repeat_count);
    }
    const ContextGenotypePriors& genotypePriors(dopt.getIndelGenotypePriors().getContextSpecificPriorSet(patternRepeatCount));

    for (const auto& score : forcedOutputIndelSampleData.read_path_lnp)
    {
        const unsigned readIndex(score.first);
        const ReadPathScores& forcedOutputAlleleReadScores(score.second);

        if (! forcedOutputAlleleReadScores.is_tier1_read) continue;

        double refAllele_lnp(forcedOutputAlleleReadScores.ref);
        const double forcedOutputAllele_lnp(forcedOutputAlleleReadScores.indel);

        double maxVariantAllele_lnp(log0);
        const unsigned maxVariantAlleleIndex(updateFromAlleleGroup(sampleId, forcedOutputIndelKey, readIndex,
                                                                   variantAlleleGroup, refAllele_lnp,
                                                                   maxVariantAllele_lnp));

        /// TEMPORARY: for now we compress REF and <*> to one state until we have a way to report this:
        if (maxVariantAllele_lnp > refAllele_lnp)
        {
            refAllele_lnp = maxVariantAllele_lnp;
        }
        maxVariantAllele_lnp = log0;

        const IndelKey* maxVariantAlleleKeyPtr(nullptr);
        if (maxVariantAlleleIndex<variantAlleleGroup.size())
        {
            maxVariantAlleleKeyPtr = &variantAlleleGroup.key(maxVariantAlleleIndex);
        }
        else
        {
            maxVariantAlleleKeyPtr = &forcedOutputIndelKey;
        }

        accumulateLogLhoodFromReadObservation(dopt, sampleOptions, forcedOutputAlleleReadScores.nsite,
                                              forcedOutputAlleleReadScores.read_length, is3AlleleModel,
                                              refAllele_lnp, forcedOutputAllele_lnp, maxVariantAllele_lnp,
                                              forcedOutputIndelKey, *maxVariantAlleleKeyPtr,
                                              logLhood);

        //---------------------------------------------------
        // get support counts for each allele
        //

        alleleLhood[0] = refAllele_lnp;
        alleleLhood[1] = forcedOutputAllele_lnp;
#if 0
        if (is3AlleleModel)
        {
            alleleLhood[2] = maxVariantAllele_lnp;
        }
#endif
        updateSupportingReadStats(dopt, readSupportTheshold, forcedOutputAlleleReadScores.nsite, forcedOutputAlleleReadScores.is_fwd_strand, alleleLhood, locusReadStats);

    }

    const bool isHaploid(groupLocusPloidy==1);
    logLhoodToLocusGenotype(genotypePriors, gtCount, is3AlleleModel, isHaploid, logLhood, locusGenotype);
}
