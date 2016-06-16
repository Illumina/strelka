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
        assert(alleleCount<=3);
    }

    /// return number of possible genotypes:
    ///
    // just hard-code in the cases we need for now...
    uint8_t
    genotypeCount() const
    {
        if (_ploidy==1)
        {
            return _alleleCount;
        }
        else if (_ploidy==2)
        {
            switch (_alleleCount)
            {
            case 1:
                return 1;
            case 2:
                return 3;
            case 3:
                return 6;
            }
        }
        assert(false && "Unexpected ploidy or haplotype count");
        return 0;
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


/// first push into a generalized Allele Genotype object which
/// could apply to SNVs and indels
///
/// to make this a high performance object, we fix the array sizes
/// for now
///
/// genotype ordering follows VCF convention
////
/// first allele is always the implicit reference allele
///
struct AlleleGroupGenotype
{
    static const uint8_t MAX_GENOTYPE_COUNT = 6;

    bool
    isNonReferenceGenotype() const
    {
        return (variantAlleleQuality != 0);
    }

    unsigned maxGenotypeIndex;
    unsigned maxGenotypeIndexPolymorphic;

    double variantAlleleQuality;
    double genotypeQuality;
    double genotypeQualityPolymorphic;
    double posteriorProb[MAX_GENOTYPE_COUNT];
    double phredLoghood[MAX_GENOTYPE_COUNT];
};



namespace AG_GENOTYPE
{
enum index_t
{
    HOMREF,
    HOM0,
    HET0,
    HOM1,
    HET1,
    HET01,
    SIZE
};

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

/// return the heterozygous genotype including alleleId and reference
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

/// return the homozygous genotype including alleleId and reference
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

}


struct GenotypePriors
{
    explicit
    GenotypePriors(
        const double theta)
    {
        static const double log0(-std::numeric_limits<double>::infinity());

        prior2AlleleDiploid[AG_GENOTYPE::HOMREF]=std::log(1.-(theta*3/2));
        prior2AlleleDiploid[AG_GENOTYPE::HOM0]=std::log(theta/2.);
        prior2AlleleDiploid[AG_GENOTYPE::HET0]=std::log(theta);

        prior2AlleleDiploidPolymorphic[AG_GENOTYPE::HOMREF]=std::log(0.25);
        prior2AlleleDiploidPolymorphic[AG_GENOTYPE::HOM0]=std::log(0.25);
        prior2AlleleDiploidPolymorphic[AG_GENOTYPE::HET0]=std::log(0.5);

        prior2AlleleHaploid[AG_GENOTYPE::HOMREF]=std::log(1.-theta);
        prior2AlleleHaploid[AG_GENOTYPE::HOM0]=std::log(theta);
        prior2AlleleHaploid[AG_GENOTYPE::HET0]=log0;

        prior2AlleleHaploidPolymorphic[AG_GENOTYPE::HOMREF]=std::log(0.5);
        prior2AlleleHaploidPolymorphic[AG_GENOTYPE::HOM0]=std::log(0.5);
        prior2AlleleHaploidPolymorphic[AG_GENOTYPE::HET0]=log0;

        prior3AlleleDiploid[AG_GENOTYPE::HOMREF]=std::log(1.-(theta*3./2.));
        prior3AlleleDiploid[AG_GENOTYPE::HOM0]=std::log(theta/2.);
        prior3AlleleDiploid[AG_GENOTYPE::HET0]=std::log(theta);
        prior3AlleleDiploid[AG_GENOTYPE::HOM1]=std::log(theta*theta/2);
        prior3AlleleDiploid[AG_GENOTYPE::HET1]=std::log(theta*theta);
        prior3AlleleDiploid[AG_GENOTYPE::HET01]=std::log(theta*theta);

        prior3AlleleDiploidPolymorphic[AG_GENOTYPE::HOMREF]=std::log(0.25);
        prior3AlleleDiploidPolymorphic[AG_GENOTYPE::HOM0]=std::log(0.25);
        prior3AlleleDiploidPolymorphic[AG_GENOTYPE::HET0]=std::log(0.5);
        prior3AlleleDiploidPolymorphic[AG_GENOTYPE::HOM1]=std::log(0.25*theta);
        prior3AlleleDiploidPolymorphic[AG_GENOTYPE::HET1]=std::log(0.5*theta);
        prior3AlleleDiploidPolymorphic[AG_GENOTYPE::HET01]=std::log(0.5*theta);

        prior3AlleleHaploid[AG_GENOTYPE::HOMREF]=std::log(1.-theta);
        prior3AlleleHaploid[AG_GENOTYPE::HOM0]=std::log(theta);
        prior3AlleleHaploid[AG_GENOTYPE::HET0]=log0;
        prior3AlleleHaploid[AG_GENOTYPE::HOM1]=std::log(theta*theta);
        prior3AlleleHaploid[AG_GENOTYPE::HET1]=log0;
        prior3AlleleHaploid[AG_GENOTYPE::HET01]=log0;

        prior3AlleleHaploidPolymorphic[AG_GENOTYPE::HOMREF]=std::log(0.5);
        prior3AlleleHaploidPolymorphic[AG_GENOTYPE::HOM0]=std::log(0.5);
        prior3AlleleHaploidPolymorphic[AG_GENOTYPE::HET0]=log0;
        prior3AlleleHaploidPolymorphic[AG_GENOTYPE::HOM1]=std::log(0.5*theta);
        prior3AlleleHaploidPolymorphic[AG_GENOTYPE::HET1]=log0;
        prior3AlleleHaploidPolymorphic[AG_GENOTYPE::HET01]=log0;

    }

    const double*
    get2Allele(const bool isHaploid) const
    {
        if (isHaploid)
        {
            return prior2AlleleHaploid;
        }
        else
        {
            return prior2AlleleDiploid;
        }
    }

    const double*
    get3Allele(const bool isHaploid) const
    {
        if (isHaploid)
        {
            return prior3AlleleHaploid;
        }
        else
        {
            return prior3AlleleDiploid;
        }
    }

    const double*
    get2AllelePolymorphic(const bool isHaploid) const
    {
        if (isHaploid)
        {
            return prior2AlleleHaploidPolymorphic;
        }
        else
        {
            return prior2AlleleDiploidPolymorphic;
        }
    }

    const double*
    get3AllelePolymorphic(const bool isHaploid) const
    {
        if (isHaploid)
        {
            return prior3AlleleHaploidPolymorphic;
        }
        else
        {
            return prior3AlleleDiploidPolymorphic;
        }
    }

    double prior2AlleleDiploid[3];
    double prior3AlleleDiploid[6];

    double prior2AlleleHaploid[3];
    double prior3AlleleHaploid[6];

    double prior2AlleleDiploidPolymorphic[3];
    double prior3AlleleDiploidPolymorphic[6];

    double prior2AlleleHaploidPolymorphic[3];
    double prior3AlleleHaploidPolymorphic[6];
};



void
getVariantAlleleGroupGenotypeLhoods(
    const starling_base_deriv_options& dopt,
    const starling_sample_options& sampleOptions,
    const GenotypePriors& genotypePriors,
    const unsigned groupLocusPloidy,
    const unsigned sampleId,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    AlleleGroupGenotype& locusGenotype);

void
getGenotypeLhoodsForForcedOutputAllele(
    const starling_base_deriv_options& dopt,
    const starling_sample_options& sampleOptions,
    const GenotypePriors& genotypePriors,
    const unsigned groupLocusPloidy,
    const unsigned sampleId,
    const OrthogonalVariantAlleleCandidateGroup& variantAlleleGroup,
    const OrthogonalVariantAlleleCandidateGroup& forcedOutputAlleleGroup,
    const unsigned forcedOutputAlleleIndex,
    AlleleGroupGenotype& locusGenotype);
