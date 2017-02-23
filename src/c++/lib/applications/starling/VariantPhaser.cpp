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
/*
 *  Author: Sangtae Kim
 */

#include "VariantPhaser.hh"

template <class T>
void VariantPhaser::processLocus(std::unique_ptr<T> locusPtr)
{
    if (not _opt.isUseCodonPhaser)
    {
        // bypass phasing
        _sink->process(std::move(locusPtr));
        return;
    }

    auto curActiveRegionId(locusPtr->getActiveRegionId());

    if (curActiveRegionId < 0)
    {
        // locus is outside of active region
        outputBuffer();
        _sink->process(std::move(locusPtr));
    }
    else if (curActiveRegionId != _activeRegionId)
    {
        // new active region
        outputBuffer();

        _activeRegionId = curActiveRegionId;
        addLocusToBuffer(std::move(locusPtr));
    }
    else
    {
        // locus is within the existing active region
        addLocusToBuffer(std::move(locusPtr));
    }
}

void
VariantPhaser::
process(std::unique_ptr<GermlineSiteLocusInfo> locusPtr)
{
    processLocus(std::move(locusPtr));
}

void
VariantPhaser::
process(std::unique_ptr<GermlineIndelLocusInfo> locusPtr)
{
    processLocus(std::move(locusPtr));
}


/// batch process all loci belonging to the same active region
void
VariantPhaser::
outputBuffer()
{
    // create phase records
    for (unsigned sampleId(0); sampleId < _sampleCount; ++sampleId)
        createPhaseRecord(sampleId);

    // move loci to the next pipe
    for (auto& locusPtr : _locusBuffer)
    {
        if (isInstanceOf<GermlineSiteLocusInfo>(*locusPtr))
        {
            _sink->process(downcast<GermlineSiteLocusInfo>(std::move(locusPtr)));
        }
        else if (isInstanceOf<GermlineIndelLocusInfo>(*locusPtr))
        {
            _sink->process(downcast<GermlineIndelLocusInfo>(std::move(locusPtr)));
        }
        else
        {
            throw std::bad_cast();
        }
    }

    _locusBuffer.clear();
}

void
VariantPhaser::
createPhaseRecord(unsigned sampleId)
{
    if (not isBuffer()) return;

    // simple hack to select 2 haplotypes assuming all variants in _locusBuffer are correct
    for (unsigned haplotypeIdPairIndex(0); haplotypeIdPairIndex<HaplotypeIdPair::SIZE; ++haplotypeIdPairIndex)
        _possibleHaplotypeIdPair[haplotypeIdPairIndex] = true;

    unsigned numHetVariants(0);
    for (const auto& locusPtr : _locusBuffer)
    {
        const auto& sampleInfo = locusPtr->getSample(sampleId);
        if (sampleInfo.isVariant())
        {
            _possibleHaplotypeIdPair[0] = false;
            const auto& maxGt = sampleInfo.max_gt();
            bool isHet = maxGt.isHet();

            // isConflict==true means that the het variant is not from the top 2 selected haplotypes
            bool isConflict(maxGt.isConflict());
            if (isHet and !isConflict)
                ++numHetVariants;

            if ((not isConflict) and (locusPtr->getActiveRegionId() >= 0))
            {
                inferHaplotypePair(maxGt.getAllele0HaplotypeId(), isHet);

                if (maxGt.getPloidy() == 2)
                    inferHaplotypePair(maxGt.getAllele1HaplotypeId(), isHet);
            }
        }
    }

    // no phasing is needed if there's 0 or 1 het variant
    if (numHetVariants <= 1)
        return;

    // select the most plausible haplotype id pair
    uint8_t selectedHaplotypeIdPairIndex(0);
    for (unsigned hapIdPairIndex(0); hapIdPairIndex<6; ++hapIdPairIndex)
    {
        if (_possibleHaplotypeIdPair[hapIdPairIndex])
        {
            selectedHaplotypeIdPairIndex = hapIdPairIndex;
            break;
        }
    }

    // to record the haplotype id of the first nonref allele of the first variant
    // i.e. to write 0|1 instead of 1|0 for the first variant of the phase block
    uint8_t haplotypeIdOfFirstNonRefAllele(0);

    for (auto& locusPtr : _locusBuffer)
    {
        auto& sampleInfo = locusPtr->getSample(sampleId);
        if (not sampleInfo.isVariant()) continue;

        sampleInfo.phaseSetId = (_activeRegionId+1);
        auto& maxGenotype = sampleInfo.max_gt();

        if ((not maxGenotype.isHet()) or maxGenotype.isConflict()) continue;

        if (selectedHaplotypeIdPairIndex != HetHap1Hap2)
        {
            // simple case where there's no flipped genotype
            sampleInfo.phaseSetId = (_activeRegionId+1);
            maxGenotype.setPhased(false);
        }
        else
        {
            auto allele0HaplotypeId(maxGenotype.getAllele0HaplotypeId());
            auto allele1HaplotypeId(maxGenotype.getAllele1HaplotypeId());

            if (haplotypeIdOfFirstNonRefAllele == 0)
            {
                if (allele0HaplotypeId > 0)
                    haplotypeIdOfFirstNonRefAllele = allele0HaplotypeId;
                else if (allele1HaplotypeId == 1)
                    haplotypeIdOfFirstNonRefAllele = 2;
                else if (allele1HaplotypeId == 2)
                    haplotypeIdOfFirstNonRefAllele = 1;
            }

            bool isFlip(false);
            if (allele0HaplotypeId == 0)
            {
                if (allele1HaplotypeId == haplotypeIdOfFirstNonRefAllele)
                    isFlip = true;
            }
            else
            {
                if (allele0HaplotypeId != haplotypeIdOfFirstNonRefAllele)
                    isFlip = true;
            }
            maxGenotype.setPhased(isFlip);
        }
    }
}


/// This method is called for each non-conflict variant
/// \param haplotypeId 0: no variant, 1: variant is from hap1, 2: variant is from hap2, 3: variant is from both hap1 and hap2
/// \param isHet true if the variant is heterozygous
void
VariantPhaser::inferHaplotypePair(const uint8_t haplotypeId, const bool isHet)
{
    if (haplotypeId == 0) return;


    switch (haplotypeId)
    {
        case 1:
            // the variant is lying on haplotype 1
            if (isHet)
            {
                // impossible haplotype id pairs if the variant is correct
                _possibleHaplotypeIdPair[HetHap2] = false;
                _possibleHaplotypeIdPair[HomHap1] = false;
                _possibleHaplotypeIdPair[HomHap2] = false;
            }
            else
            {
                _possibleHaplotypeIdPair[HetHap1] = false;
                _possibleHaplotypeIdPair[HetHap2] = false;
                _possibleHaplotypeIdPair[HetHap1Hap2] = false;
                _possibleHaplotypeIdPair[HomHap2] = false;
            }
            break;
        case 2:
            // the variant is lying on haplotype 2
            if (isHet)
            {
                _possibleHaplotypeIdPair[HetHap1] = false;
                _possibleHaplotypeIdPair[HomHap1] = false;
                _possibleHaplotypeIdPair[HomHap2] = false;
            }
            else
            {
                _possibleHaplotypeIdPair[HetHap1] = false;
                _possibleHaplotypeIdPair[HetHap2] = false;
                _possibleHaplotypeIdPair[HetHap1Hap2] = false;
                _possibleHaplotypeIdPair[HomHap1] = false;
            }
            break;
        case 3:
            // the variant is lying on both haplotype 1 and 2
            if (isHet)
            {
                _possibleHaplotypeIdPair[HetHap1Hap2] = false;
                _possibleHaplotypeIdPair[HomHap1] = false;
                _possibleHaplotypeIdPair[HomHap2] = false;
            }
            else
            {
                _possibleHaplotypeIdPair[HetHap1] = false;
                _possibleHaplotypeIdPair[HetHap2] = false;
            }
            break;
        case 4:
            // do nothing
            break;
        default:
            assert(false);
    }
}

void
VariantPhaser::
flush_impl()
{
    if (_opt.isUseCodonPhaser)
    {
        outputBuffer();
    }
}
