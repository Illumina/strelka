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
 *
 *  Created on: Sep 10, 2013
 *  Author: Morten Kallberg
 */

#include "codon_phaser.hh"

#include <array>
#include <functional>
#include <vector>
#include <iostream>

//#define DEBUG_CODON



#ifdef DEBUG_CODON
#include "blt_util/log.hh"
#endif

template <class T>
void Codon_phaser::processLocus(std::unique_ptr<T> locusPtr)
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
Codon_phaser::
process(std::unique_ptr<GermlineSiteLocusInfo> locusPtr)
{
    processLocus(std::move(locusPtr));
}

void
Codon_phaser::
process(std::unique_ptr<GermlineIndelLocusInfo> locusPtr)
{
    processLocus(std::move(locusPtr));
}


void
Codon_phaser::
outputBuffer()
{
    for (unsigned sampleId(0); sampleId < _sampleCount; ++sampleId)
        createPhasedRecord(sampleId);

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

static void inferComplexIds(
        const uint8_t complexAlleleId,
        const bool isHet,
        std::vector<bool>& possibleComplexAlleleIdPair)
{
    if (complexAlleleId == 0) return;

    // select 2 complex allele IDs: (0,1)=>1, (0,2)=>2, (1,2)=>3, (1,1)=>4, or (2,2)=>5
    switch (complexAlleleId)
    {
        case 1:
            // haplotype 1 only
            if (isHet)
            {
                possibleComplexAlleleIdPair[2] = false;
                possibleComplexAlleleIdPair[4] = false;
                possibleComplexAlleleIdPair[5] = false;
            }
            else
            {
                possibleComplexAlleleIdPair[1] = false;
                possibleComplexAlleleIdPair[2] = false;
                possibleComplexAlleleIdPair[3] = false;
                possibleComplexAlleleIdPair[5] = false;
            }
            break;
        case 2:
            // haplotype 2 only
            if (isHet)
            {
                possibleComplexAlleleIdPair[1] = false;
                possibleComplexAlleleIdPair[4] = false;
                possibleComplexAlleleIdPair[5] = false;
            }
            else
            {
                possibleComplexAlleleIdPair[1] = false;
                possibleComplexAlleleIdPair[2] = false;
                possibleComplexAlleleIdPair[3] = false;
                possibleComplexAlleleIdPair[4] = false;
            }
            break;
        case 3:
            // both haplotype 1 and 2
            if (isHet)
            {
                possibleComplexAlleleIdPair[3] = false;
                possibleComplexAlleleIdPair[4] = false;
                possibleComplexAlleleIdPair[5] = false;
            }
            else
            {
                possibleComplexAlleleIdPair[1] = false;
                possibleComplexAlleleIdPair[2] = false;
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
Codon_phaser::
createPhasedRecord(unsigned sampleId)
{
    if (not isBuffer()) return;

    // select 2 complex allele IDs: (0,1)=>1, (0,2)=>2, (1,2)=>3, (1,1)=>4, or (2,2)=>5

    std::vector<bool> possibleComplexAlleleIdPair(6);
    for (unsigned complexAlleleIdPairIndex(0); complexAlleleIdPairIndex<6; ++complexAlleleIdPairIndex)
        possibleComplexAlleleIdPair[complexAlleleIdPairIndex] = true;

    unsigned numHetVariants(0);
    for (const auto& locusPtr : _locusBuffer)
    {
        const auto& sampleInfo = locusPtr->getSample(sampleId);
        if (sampleInfo.isVariant())
        {
            possibleComplexAlleleIdPair[0] = false;
            const auto& maxGt = sampleInfo.max_gt();
            bool isHet = maxGt.isHet();
            bool isConflict(maxGt.isConflict());
            if (isHet and !isConflict)
                ++numHetVariants;

            if (isHet and (not isConflict) and (locusPtr->getActiveRegionId() >= 0))
            {
                inferComplexIds(maxGt.getAllele0ComplexAlleleId(), isHet, possibleComplexAlleleIdPair);

                if (maxGt.getPloidy() == 2)
                    inferComplexIds(maxGt.getAllele1ComplexAlleleId(), isHet, possibleComplexAlleleIdPair);
            }
        }
    }

    if (numHetVariants <= 1)
        return;

        // select 2 complex alleles: (0,1)=>1, (0,2)=>2, (1,2)=>3, (1,1)=>4, or (2,2)=>5
        uint8_t selectedComplexAllelePairIndex(0);
        if (not possibleComplexAlleleIdPair[0])
        {
            for (unsigned complexAlleleIdPairIndex(1); complexAlleleIdPairIndex<6; ++complexAlleleIdPairIndex)
            {
                if (possibleComplexAlleleIdPair[complexAlleleIdPairIndex])
                {
                    selectedComplexAllelePairIndex = complexAlleleIdPairIndex;
                    break;
                }
            }
        }

    uint8_t firstComplexAlleleId(0);

    for (auto& locusPtr : _locusBuffer)
    {
        auto& sampleInfo = locusPtr->getSample(sampleId);
        if (sampleInfo.isVariant())
        {
            sampleInfo.phaseSetId = (_activeRegionId+1);
            auto& maxGenotype = sampleInfo.max_gt();
            auto allele0ComplexAlleleId(maxGenotype.getAllele0ComplexAlleleId());
            if (maxGenotype.isHet() and (not maxGenotype.isConflict()) and (maxGenotype.getAllele1ComplexAlleleId() != 4))
            {
                auto allele1ComplexAlleleId(maxGenotype.getAllele1ComplexAlleleId());

                if (selectedComplexAllelePairIndex == 3 and firstComplexAlleleId == 0)
                {
                    if (allele0ComplexAlleleId > 0)
                        firstComplexAlleleId = allele0ComplexAlleleId;
                    else if (allele1ComplexAlleleId == 1)
                        firstComplexAlleleId = 2;
                    else if (allele1ComplexAlleleId == 2)
                        firstComplexAlleleId = 1;
                }

                bool isFlip(false);
                if (selectedComplexAllelePairIndex == 3)
                {
                    if (allele0ComplexAlleleId == 0)
                    {
                        if (allele1ComplexAlleleId == firstComplexAlleleId)
                            isFlip = true;
                    }
                    else
                    {
                        if (allele0ComplexAlleleId != firstComplexAlleleId)
                            isFlip = true;
                    }
                }
                maxGenotype.setPhased(isFlip);
            }
        }
    }
}

void
Codon_phaser::
flush_impl()
{
    if (_opt.isUseCodonPhaser)
    {
        outputBuffer();
    }
}
