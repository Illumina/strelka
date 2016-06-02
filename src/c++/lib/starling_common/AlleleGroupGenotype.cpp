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



void
getVariantAlleleGroupGenotypeLhoods(
    const starling_base_deriv_options &dopt,
    const GenotypePriors& genotypePriors,
    const unsigned sampleId,
    const OrthogonalVariantAlleleCandidateGroup &alleleGroup,
    AlleleGroupGenotype &locusGenotype)
{
    // fix this on first pass
    static const uint8_t ploidy(2);

    assert(ploidy>0);
    assert(alleleGroup.size()<=ploidy);
    assert(ploidy<3);

    const uint8_t nonRefAlleleCount(alleleGroup.size());
    const uint8_t fullAlleleCount(nonRefAlleleCount+1);

    GenotypeInfo ginfo(ploidy,fullAlleleCount);
    const unsigned gtCount(ginfo.genotypeCount());
    assert(gtCount<=AlleleGroupGenotype::MAX_GENOTYPE_COUNT);
    double logLhood[AlleleGroupGenotype::MAX_GENOTYPE_COUNT];
    std::fill(logLhood,logLhood+gtCount,0.);

    const unsigned nonRefAllele0Index(0);
//    const indel_key& ik0(hg.key(nonRefAlleleIndex0));
    const IndelData& id0(alleleGroup.data(nonRefAllele0Index));
    const IndelSampleData& isd0(id0.getSampleData(sampleId));

    const bool is3AlleleModel(fullAlleleCount==3);
    const unsigned nonRefAllele1Index(is3AlleleModel ? 1 : 0);
    const IndelData& id1(alleleGroup.data(nonRefAllele1Index));
    const IndelSampleData& isd1(id1.getSampleData(sampleId));

    for (const auto& score : isd0.read_path_lnp)
    {
        const auto readIndex(score.first);
        const ReadPathScores& Allele0ReadScores(score.second);

        if (! Allele0ReadScores.is_tier1_read) continue;

        const ReadPathScores* Allele1ReadScoresPtr(nullptr);
        if (is3AlleleModel)
        {
            /// skip this read if it's not aligned to both indels:
            const auto Allele1ReadScoresIter(isd1.read_path_lnp.find(readIndex));
            if (Allele1ReadScoresIter == isd1.read_path_lnp.end()) continue;
            Allele1ReadScoresPtr = (&(Allele1ReadScoresIter->second));
        }

        double refAllele_lnp(Allele0ReadScores.ref);
        const double variantAllele0_lnp(Allele0ReadScores.indel);
        double variantAllele1_lnp(variantAllele0_lnp);

        if (is3AlleleModel)
        {
            const ReadPathScores& Allele1ReadScores(*Allele1ReadScoresPtr);
            refAllele_lnp = std::max(refAllele_lnp,(double)Allele1ReadScores.ref);
            variantAllele1_lnp = Allele1ReadScores.indel;
        }

        logLhood[AG_GENOTYPE::HOMREF] += integrate_out_sites(dopt,Allele0ReadScores.nsite,refAllele_lnp,false);
        logLhood[AG_GENOTYPE::HOM0] += integrate_out_sites(dopt,Allele0ReadScores.nsite,variantAllele0_lnp,false);

        static const double loghalf(std::log(0.5));

        /// disable indel specific ratios for first pass at this (and maybe forever???)
#if 0
        double log_ref_prob(loghalf);
        double log_indel_prob(loghalf);
        if (! is_breakpoint)
        {
            static const double het_allele_ratio(0.5);
            get_het_observed_allele_ratio(path_lnp.read_length,sample_opt.min_read_bp_flank,
                                          ik,het_allele_ratio,log_ref_prob,log_indel_prob);
        }
        const double het_lnp(log_sum(refAllele_lnp+log_ref_prob,hom_lnp+log_indel_prob));
#endif

        const double het0_lnp(log_sum(refAllele_lnp+loghalf,variantAllele0_lnp+loghalf));
        logLhood[AG_GENOTYPE::HET0] += integrate_out_sites(dopt,Allele0ReadScores.nsite,het0_lnp,false);

        if (is3AlleleModel)
        {
            logLhood[AG_GENOTYPE::HOM1] += integrate_out_sites(dopt,Allele0ReadScores.nsite,variantAllele1_lnp,false);

            const double het1_lnp(log_sum(refAllele_lnp+loghalf,variantAllele1_lnp+loghalf));
            logLhood[AG_GENOTYPE::HET1] += integrate_out_sites(dopt,Allele0ReadScores.nsite,het1_lnp,false);

            const double het01_lnp(log_sum(variantAllele0_lnp+loghalf,variantAllele1_lnp+loghalf));
            logLhood[AG_GENOTYPE::HET01] += integrate_out_sites(dopt,Allele0ReadScores.nsite,het01_lnp,false);
        }
    }

    // mult by prior distro to get unnormalized pprob:
    const double* genotypeLogPrior(genotypePriors.prior2AlleleDiploid);
    if (is3AlleleModel)
    {
        genotypeLogPrior = (genotypePriors.prior3AlleleDiploid);
    }

    for (unsigned gt(0); gt<gtCount; ++gt)
    {
        locusGenotype.posteriorProb[gt] = logLhood[gt] + genotypeLogPrior[gt];
    }

    normalize_ln_distro(std::begin(locusGenotype.posteriorProb),std::begin(locusGenotype.posteriorProb)+gtCount,locusGenotype.maxGenotypeIndex);

    locusGenotype.variantAlleleQuality=error_prob_to_qphred(locusGenotype.posteriorProb[AG_GENOTYPE::HOMREF]);
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
        const double* genotypeLogPriorPolymorphic(genotypePriors.prior2AlleleDiploidPolymorphic);
        if (is3AlleleModel)
        {
            genotypeLogPriorPolymorphic = (genotypePriors.prior3AlleleDiploidPolymorphic);
        }

        for (unsigned gt(0); gt<gtCount; ++gt)
        {
            logLhood[gt] += genotypeLogPriorPolymorphic[gt];
        }
    }
    normalize_ln_distro(logLhood,logLhood+gtCount,locusGenotype.maxGenotypeIndexPolymorphic);
    locusGenotype.genotypeQualityPolymorphic=error_prob_to_qphred(prob_comp(logLhood,logLhood+gtCount,locusGenotype.maxGenotypeIndexPolymorphic));
}
