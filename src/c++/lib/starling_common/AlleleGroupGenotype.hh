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
