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
/*
 *  Author: Sangtae Kim
 */

#include "VariantPhaser.hh"

template <class T>
void VariantPhaser::processLocus(std::unique_ptr<T> locusPtr)
{
    if (not _opt.isUseVariantPhaser)
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

    // simple hack for checking whether 2 alternative haplotypes (haplotypeId 1 and 2) are selected
    auto isHetHap1(false);
    auto isHetHap2(false);
    unsigned numHetVariants(0);
    for (const auto& locusPtr : _locusBuffer)
    {
        const auto& sampleInfo = locusPtr->getSample(sampleId);
        if (sampleInfo.isVariant())
        {
            const auto& maxGenotype = sampleInfo.max_gt();
            bool isHet = maxGenotype.isHet();

            // isConflict==true means that the het variant is not from the top 2 selected haplotypes
            bool isConflict(maxGenotype.isConflict());
            if (isHet and !isConflict)
                ++numHetVariants;

            if (isHet and (not isConflict) and (locusPtr->getActiveRegionId() >= 0))
            {
                if (maxGenotype.getAllele0HaplotypeId() == 1)
                    isHetHap1 = true;
                else if (maxGenotype.getAllele0HaplotypeId() == 2)
                    isHetHap2 = true;

                if (maxGenotype.getPloidy() == 2)
                {
                    if (maxGenotype.getAllele1HaplotypeId() == 1)
                        isHetHap1 = true;
                    else if (maxGenotype.getAllele1HaplotypeId() == 2)
                        isHetHap2 = true;
                }
            }
        }
    }

    // no phasing is needed if there's 0 or 1 het variant
    if (numHetVariants <= 1)
        return;

    // to record the haplotype id of the first nonref allele of the first variant
    // i.e. to write 0|1 instead of 1|0 for the first variant of the phase block
    uint8_t haplotypeIdOfFirstNonRefAllele(0);

    for (auto& locusPtr : _locusBuffer)
    {
        auto& sampleInfo = locusPtr->getSample(sampleId);
        if (not sampleInfo.isVariant()) continue;

        auto& maxGenotype = sampleInfo.max_gt();

        if ((not maxGenotype.isHet()) or maxGenotype.isConflict()) continue;

        sampleInfo.phaseSetId = (_activeRegionId+1);
        if ((not isHetHap1) or (not isHetHap2))
        {
            // simple case where there's no flipped genotype
            // one haplotype is the reference
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

void
VariantPhaser::
flush_impl()
{
    if (_opt.isUseVariantPhaser)
    {
        outputBuffer();
    }
}
