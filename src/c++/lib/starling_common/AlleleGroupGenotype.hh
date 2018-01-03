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

#pragma once

#include "starling_base_shared.hh"
#include "starling_common/OrthogonalVariantAlleleCandidateGroup.hh"
#include "LocusSupportingReadStats.hh"

#include "htsapi/vcf_util.hh"


namespace AG_GENOTYPE
{
/// genotype enum generalized to handle-multi sample as follows:
enum index_t
{
    HOMREF,
    HET0, ///< ref + most likely alt
    HOM0, ///< most likely alt
    HET1, ///< ref + any other alt
    HOM1, ///< any other alt
    HET01, ///< most likely alt + any other alt
    SIZE
};
}


/// Useful genotype priors generated for a given theta value
struct ContextGenotypePriors
{
    void
    initialize(
        const double theta)
    {
        static const double log0(-std::numeric_limits<double>::infinity());
        /// Penalty for calling a genotype that includes allele 1 but excludes allele 0 (which has more support)
        static const double allele0SkipPenalty(theta);
        // TODO: (1) add 1-allele0SkipPenalty factor when allele 0 is present and allele 1 is absent
        //       (2) adjust homref probabilities such that everything sums to 1
        //       (3) experiment with the value of allele0SkipPenalty

        _isInitialized = true;

        priorNAlleleDiploid[AG_GENOTYPE::HOMREF] = std::log(1. - (theta * 3. / 2.));
        priorNAlleleDiploid[AG_GENOTYPE::HOM0] = std::log(theta / 2.);
        priorNAlleleDiploid[AG_GENOTYPE::HET0] = std::log(theta);
        priorNAlleleDiploid[AG_GENOTYPE::HOM1] = std::log(theta * allele0SkipPenalty / 2);
        priorNAlleleDiploid[AG_GENOTYPE::HET1] = std::log(theta * allele0SkipPenalty);
        priorNAlleleDiploid[AG_GENOTYPE::HET01] = std::log(theta * theta);

        priorNAlleleDiploidPolymorphic[AG_GENOTYPE::HOMREF] = std::log(0.25);
        priorNAlleleDiploidPolymorphic[AG_GENOTYPE::HOM0] = std::log(0.25);
        priorNAlleleDiploidPolymorphic[AG_GENOTYPE::HET0] = std::log(0.5);
        priorNAlleleDiploidPolymorphic[AG_GENOTYPE::HOM1] = std::log(0.25 * allele0SkipPenalty);
        priorNAlleleDiploidPolymorphic[AG_GENOTYPE::HET1] = std::log(0.5 * allele0SkipPenalty);
        priorNAlleleDiploidPolymorphic[AG_GENOTYPE::HET01] = std::log(0.5 * theta);

        priorNAlleleHaploid[AG_GENOTYPE::HOMREF] = std::log(1. - theta);
        priorNAlleleHaploid[AG_GENOTYPE::HOM0] = std::log(theta);
        priorNAlleleHaploid[AG_GENOTYPE::HET0] = log0;
        priorNAlleleHaploid[AG_GENOTYPE::HOM1] = std::log(theta * allele0SkipPenalty);
        priorNAlleleHaploid[AG_GENOTYPE::HET1] = log0;
        priorNAlleleHaploid[AG_GENOTYPE::HET01] = log0;

        priorNAlleleHaploidPolymorphic[AG_GENOTYPE::HOMREF] = std::log(0.5);
        priorNAlleleHaploidPolymorphic[AG_GENOTYPE::HOM0] = std::log(0.5);
        priorNAlleleHaploidPolymorphic[AG_GENOTYPE::HET0] = log0;
        priorNAlleleHaploidPolymorphic[AG_GENOTYPE::HOM1] = std::log(0.5 * allele0SkipPenalty);
        priorNAlleleHaploidPolymorphic[AG_GENOTYPE::HET1] = log0;
        priorNAlleleHaploidPolymorphic[AG_GENOTYPE::HET01] = log0;
    }

    const double*
    getNAllele(const bool isHaploid) const
    {
        assert(_isInitialized);

        if (isHaploid)
        {
            return priorNAlleleHaploid;
        }
        else
        {
            return priorNAlleleDiploid;
        }
    }

    const double*
    getNAllelePolymorphic(const bool isHaploid) const
    {
        assert(_isInitialized);

        if (isHaploid)
        {
            return priorNAlleleHaploidPolymorphic;
        }
        else
        {
            return priorNAlleleDiploidPolymorphic;
        }
    }

    double priorNAlleleDiploid[AG_GENOTYPE::SIZE];
    double priorNAlleleHaploid[AG_GENOTYPE::SIZE];
    double priorNAlleleDiploidPolymorphic[AG_GENOTYPE::SIZE];
    double priorNAlleleHaploidPolymorphic[AG_GENOTYPE::SIZE];
private:
    bool _isInitialized = false;
};



struct GenotypePriorSet
{

    explicit GenotypePriorSet(
        const std::string& thetaJsonFilename);

    std::map<unsigned, std::vector<double> >
    initializeThetas();

    void
    initializePriors(
        const std::map<unsigned, std::vector<double> >& thetas);

    const ContextGenotypePriors&
    getContextSpecificPriorSet(
        unsigned repeatingPatternSize,
        const unsigned patternRepeatCount) const
    {
        assert(repeatingPatternSize>0);
        assert(patternRepeatCount>0);

        // STR pattern lengths that we don't represent are treated as the homopolymer length
        if (repeatingPatternSize > _priors.size())
        {
            repeatingPatternSize = 1;
        }

        const unsigned repeatingPatternSizeIndex(repeatingPatternSize-1);

        // STR repeat counts we don't represent are treated as the highest repeat value that we do represent:
        const auto& STRPatternPriors(_priors[repeatingPatternSizeIndex]);
        const unsigned patternRepeatCountIndex(std::min(patternRepeatCount,static_cast<unsigned>(STRPatternPriors.size()))-1);
        return STRPatternPriors[patternRepeatCountIndex];
    }

private:
    std::vector<std::vector<ContextGenotypePriors>> _priors;
};


/// contrast group contains alleles intended for an "other" category, such as reported by the <*> allele in
/// vcf
void
getVariantAlleleGroupGenotypeLhoodsForSample(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const starling_sample_options& sampleOptions,
    const unsigned callerPloidy,
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const OrthogonalVariantAlleleCandidateGroup& contrastGroup,
    std::vector<double>& genotypeLogLhood,
    LocusSupportingReadStats& locusReadStats);
