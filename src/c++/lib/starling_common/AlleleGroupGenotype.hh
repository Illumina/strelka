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

#pragma once

#include "starling_base_shared.hh"
#include "starling_common/OrthogonalVariantAlleleCandidateGroup.hh"
#include "LocusSupportingReadStats.hh"

#include "htsapi/vcf_util.hh"


/// produce various genotype information based on
/// - ploidy
/// - no of haplotypes
///
struct GenotypeInfo
{
    GenotypeInfo(
        const uint8_t ploidy,
        const uint8_t alleleCount)
        : _ploidy(ploidy), _alleleCount(alleleCount)
    {
        assert(ploidy>0);
        assert(ploidy<=2);
        assert(alleleCount>0);
    }

    /// return number of possible genotypes:
    ///
    // just hard-code in the cases we need for now...
    uint8_t
    genotypeCount() const
    {
        return VcfGenotypeUtil::getGenotypeCount(_ploidy, _alleleCount);
    }

#if 0
    /// for a specific genotype number, the haplotype
    /// for chomosomeIndex N, where N<ploidy
    uint8_t
    getHapIndexFromGenotypeIndex(
        const uint8_t genotypeIndex,
        const uint8_t chromIndex) const
    {

    }

    uint8_t
    getIndexFromAlleles(
        const uint8_t a1,
        const uint8_t a2) const
    {
        assert(a1<=a2);
        return (a2*(a2+1)/2)+a1;
    }
#endif

private:
    uint8_t _ploidy;
    uint8_t _alleleCount;
};


/// generalized allele genotype object which could apply to SNVs and indels
///
/// to make this a high performance object, we fix the array sizes
/// for now
///
/// genotype ordering follows VCF convention
///
/// first allele is always the implicit reference allele
///
struct AlleleGroupGenotype
{
    static const uint8_t MAX_GENOTYPE_COUNT = 6;

    bool
    isNonReferenceGenotype() const
    {
        return (anyVariantAlleleQuality != 0);
    }

    unsigned maxGenotypeIndex;
    unsigned maxGenotypeIndexPolymorphic;

    double anyVariantAlleleQuality;
    double genotypeQuality;
    double genotypeQualityPolymorphic;
    double posteriorProb[MAX_GENOTYPE_COUNT];
    double phredLoghood[MAX_GENOTYPE_COUNT];
};



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

#if 0
/// ploidy is infered from command-line arguments:
inline
index_t
getGenotypeId(
    const unsigned alleleId0)
{
    assert(alleleId0 <= 2);

    if (alleleId0 == 0)
    {
        return HOMREF;
    }
    else if (alleleId0 == 1)
    {
        return HOM0;
    }
    else if (alleleId0 == 2)
    {
        return HOM1;
    }
    assert(false and "Shouldn't be here");
    return HOMREF;
}

inline
index_t
getGenotypeId(
    const unsigned alleleId0,
    const unsigned alleleId1)
{
    assert(alleleId0 <= 2);
    assert(alleleId1 <= 2);
    assert(alleleId0 <= alleleId1);

    if (alleleId0 == 0)
    {
        if (alleleId1 == 0)
        {
            return HOMREF;
        }
        else if (alleleId1 == 1)
        {
            return HET0;
        }
        else if (alleleId1 == 2)
        {
            return HET1;
        }
    }
    else if (alleleId0 == 1)
    {
        if (alleleId1 == 1)
        {
            return HOM0;
        }
        else if (alleleId1 == 2)
        {
            return HET01;
        }
    }
    else if (alleleId0 == 2)
    {
        if (alleleId1 == 2)
        {
            return HOM1;
        }
    }
    assert(false and "Shouldn't be here");
    return HOMREF;
}

inline
bool
isAllelePresent(
    const unsigned genotypeId,
    const unsigned alleleId)
{
    if (alleleId==0)
    {
        switch (static_cast<index_t>(genotypeId))
        {
        case HOM0:
        case HET0:
        case HET01:
            return true;
        default:
            return false;
        }
    }
    else if (alleleId == 1)
    {
        switch (static_cast<index_t>(genotypeId))
        {
        case HOM1:
        case HET1:
        case HET01:
            return true;
        default:
            return false;
        }
    }
    else
    {
        assert(false and "Unsupported alleleId");
        return false;
    }
}

/// return the heterozygous genotype composed of (alleleId, reference allele)
inline
index_t
getAlleleHetId(
    const unsigned alleleId)
{
    if (alleleId == 0)
    {
        return HET0;
    }
    else if (alleleId == 1)
    {
        return HET1;
    }
    else
    {
        assert(false and "Unsupported alleleId");
        return SIZE;
    }
}

/// return the homozygous genotype composed by alleleId
inline
index_t
getAlleleHomId(
    const unsigned alleleId)
{
    if (alleleId == 0)
    {
        return HOM0;
    }
    else if (alleleId == 1)
    {
        return HOM1;
    }
    else
    {
        assert(false and "Unsupported alleleId");
        return SIZE;
    }
}
#endif
}



struct ContextGenotypePriors
{
    void
    initialize(
        const double theta)
    {
        static const double log0(-std::numeric_limits<double>::infinity());

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
};



struct GenotypePriorSet
{
    GenotypePriorSet(
        const double lowRepeatTheta,
        const double highRepeatTheta,
        const unsigned highRepeatCount)
        : _priors(highRepeatCount)
    {
        assert(highRepeatCount>0);

        const unsigned highRepeatCountIndex(highRepeatCount-1);
        for (unsigned patternRepeatCount=1; patternRepeatCount <= highRepeatCount; ++patternRepeatCount)
        {
            const unsigned patternRepeatCountIndex(patternRepeatCount-1);
            const double highValueFraction(std::min(patternRepeatCountIndex,highRepeatCountIndex)/static_cast<double>(highRepeatCountIndex));
            const double theta((1.-highValueFraction)*lowRepeatTheta + highValueFraction*highRepeatTheta);
            _priors[patternRepeatCountIndex].initialize(theta);
        }
    }

    const ContextGenotypePriors&
    getContextSpecificPriorSet(
        const unsigned patternRepeatCount) const
    {
        assert(patternRepeatCount>0);
        const unsigned patternRepeatCountIndex(std::min(patternRepeatCount,static_cast<unsigned>(_priors.size()))-1);
        return _priors[patternRepeatCountIndex];
    }

private:
    std::vector<ContextGenotypePriors> _priors;
};


void
getVariantAlleleGroupGenotypeLhoodsForSample(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const starling_sample_options& sampleOptions,
    const unsigned callerPloidy,
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    std::vector<double>& genotypeLogLhood,
    LocusSupportingReadStats& locusReadStats);

void
getGenotypeLhoodsForForcedOutputAllele(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const starling_sample_options& sampleOptions,
    const unsigned groupLocusPloidy,
    const unsigned sampleId,
    const OrthogonalVariantAlleleCandidateGroup& variantAlleleGroup,
    const OrthogonalVariantAlleleCandidateGroup& forcedOutputAlleleGroup,
    const unsigned forcedOutputAlleleIndex,
    AlleleGroupGenotype& locusGenotype,
    LocusSupportingReadStats& locusReadStats);
