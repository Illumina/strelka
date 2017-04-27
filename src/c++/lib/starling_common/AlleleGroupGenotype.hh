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

        _isInitialized = true;

        priorNAlleleDiploid[AG_GENOTYPE::HOMREF] = std::log(1. - (theta * 3. / 2.));
        priorNAlleleDiploid[AG_GENOTYPE::HOM0] = std::log(theta / 2.);
        priorNAlleleDiploid[AG_GENOTYPE::HET0] = std::log(theta);
        priorNAlleleDiploid[AG_GENOTYPE::HOM1] = std::log(theta * theta / 2);
        priorNAlleleDiploid[AG_GENOTYPE::HET1] = std::log(theta * theta);
        priorNAlleleDiploid[AG_GENOTYPE::HET01] = std::log(theta * theta);

        priorNAlleleDiploidPolymorphic[AG_GENOTYPE::HOMREF] = std::log(0.25);
        priorNAlleleDiploidPolymorphic[AG_GENOTYPE::HOM0] = std::log(0.25);
        priorNAlleleDiploidPolymorphic[AG_GENOTYPE::HET0] = std::log(0.5);
        priorNAlleleDiploidPolymorphic[AG_GENOTYPE::HOM1] = std::log(0.25 * theta);
        priorNAlleleDiploidPolymorphic[AG_GENOTYPE::HET1] = std::log(0.5 * theta);
        priorNAlleleDiploidPolymorphic[AG_GENOTYPE::HET01] = std::log(0.5 * theta);

        priorNAlleleHaploid[AG_GENOTYPE::HOMREF] = std::log(1. - theta);
        priorNAlleleHaploid[AG_GENOTYPE::HOM0] = std::log(theta);
        priorNAlleleHaploid[AG_GENOTYPE::HET0] = log0;
        priorNAlleleHaploid[AG_GENOTYPE::HOM1] = std::log(theta * theta);
        priorNAlleleHaploid[AG_GENOTYPE::HET1] = log0;
        priorNAlleleHaploid[AG_GENOTYPE::HET01] = log0;

        priorNAlleleHaploidPolymorphic[AG_GENOTYPE::HOMREF] = std::log(0.5);
        priorNAlleleHaploidPolymorphic[AG_GENOTYPE::HOM0] = std::log(0.5);
        priorNAlleleHaploidPolymorphic[AG_GENOTYPE::HET0] = log0;
        priorNAlleleHaploidPolymorphic[AG_GENOTYPE::HOM1] = std::log(0.5 * theta);
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
    GenotypePriorSet(
        const double /*lowRepeatTheta*/,
        const double /*highRepeatTheta*/,
        const unsigned /*highRepeatCount*/)
    {
        static const unsigned highHpolRepeatCount(16);
        static const double hpolTheta[] =
        {
            0.000120268,
            5.97777E-05,
            0.000124648,
            0.000260759,
            0.000589544,
            0.002394583,
            0.007417864,
            0.022660355,
            0.04670561,
            0.082031233,
            0.124548518,
            0.149765438,
            0.168051826,
            0.187346626,
            0.207339703,
            0.225843098,
            0.248849306,
            0.27106361,
            0.334718891,
            0.348811678
        };

        static const unsigned hpolThetaSize = sizeof(hpolTheta)/sizeof(double);
        assert(hpolThetaSize >= highHpolRepeatCount);

        static const unsigned highDinucRepeatCount(9);
        static const double dinucTheta[] =
        {
            0.000120268,
            8.73757E-05,
            0.000479319,
            0.002678401,
            0.012194565,
            0.03162284,
            0.060846617,
            0.108263861,
            0.163510548,
            0.204456064,
            0.23462438,
            0.267919304,
            0.290588942,
            0.355588567,
            0.369478351,
            0.378290471,
            0.38555006,
            0.393439865,
            0.395844077,
            0.4
        };

        static const unsigned dinucThetaSize = sizeof(dinucTheta)/sizeof(double);
        assert(dinucThetaSize >= highDinucRepeatCount);

        static const unsigned maxRepeatingPatternSize(2);

        _priors.resize(maxRepeatingPatternSize);
        for (unsigned repeatingPatternSize(1); repeatingPatternSize <= maxRepeatingPatternSize; ++repeatingPatternSize)
        {
            const unsigned repeatingPatternSizeIndex(repeatingPatternSize-1);
            auto& strPatternPriors(_priors[repeatingPatternSizeIndex]);

            if (repeatingPatternSize == 1)
            {
                strPatternPriors.resize(highHpolRepeatCount);
                for (unsigned patternRepeatCount(1); patternRepeatCount <= highHpolRepeatCount; ++patternRepeatCount)
                {
                    const unsigned patternRepeatCountIndex(patternRepeatCount - 1);
                    const double theta(hpolTheta[patternRepeatCountIndex]);
                    strPatternPriors[patternRepeatCountIndex].initialize(theta);
                }
            }
            else if (repeatingPatternSize == 2)
            {
                strPatternPriors.resize(highDinucRepeatCount);
                for (unsigned patternRepeatCount(1); patternRepeatCount <= highDinucRepeatCount; ++patternRepeatCount)
                {
                    const unsigned patternRepeatCountIndex(patternRepeatCount - 1);
                    const double theta(dinucTheta[patternRepeatCountIndex]);
                    strPatternPriors[patternRepeatCountIndex].initialize(theta);
                }
            }
        }
    }

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
