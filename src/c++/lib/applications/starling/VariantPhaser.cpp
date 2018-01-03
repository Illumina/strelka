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

    // The locus is moved to the buffer unless it is outside of an active region
    auto curActiveRegionId(locusPtr->getActiveRegionId());

    if (curActiveRegionId < 0)
    {
        // locus is outside of active region in all samples
        // Create phase records for the loci in the buffer
        outputBuffer();
        _sink->process(std::move(locusPtr));
    }
    else
    {
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
    for (unsigned sampleIndex(0); sampleIndex < _sampleCount; ++sampleIndex)
        createPhaseRecord(sampleIndex);

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
createPhaseRecord(unsigned sampleIndex)
{
    // Create phase records for the specified sample.
    // Note that the locus buffer may contain loci from different active regions
    // even in the same sample
    // E.g.  ------------- AR1 -------------      (Sample 1)
    //           ---- AR2 ---   ----- AR3 ----    (Sample 2)
    //       ---------------------------------    (Loci in the buffer)
    if (not isBuffer()) return;

    // simple hack for checking whether 2 alternative haplotypes (haplotypeId 1 and 2) are selected
    // Use vectors to separate loci from different ARs
    std::vector<bool> isHetHap1;
    std::vector<bool> isHetHap2;
    std::vector<unsigned> numHetVariants;

    ActiveRegionId activeRegionId = -1;
    for (const auto& locusPtr : _locusBuffer)
    {
        const auto& sampleInfo = locusPtr->getSample(sampleIndex);
        const auto curActiveRegionId(sampleInfo.getActiveRegionId());

        if ((!sampleInfo.isVariant()) || (curActiveRegionId < 0)) continue;

        if (curActiveRegionId != activeRegionId)
        {
            // new active region; clear the vectors
            isHetHap1.push_back(false);
            isHetHap2.push_back(false);
            numHetVariants.push_back(0);

            activeRegionId = curActiveRegionId;
        }

        const auto& maxGenotype = sampleInfo.max_gt();
        const bool isHet = maxGenotype.isHet();

        // isConflict==true means that the het variant is not from the top 2 selected haplotypes
        const bool isConflict(maxGenotype.isConflict());

        if (isHet and (!isConflict))
        {
            ++(numHetVariants.back());
            if (maxGenotype.getAllele0HaplotypeId() == 1)
                isHetHap1.back() = true;
            else if (maxGenotype.getAllele0HaplotypeId() == 2)
                isHetHap2.back() = true;

            if (maxGenotype.getPloidy() == 2)
            {
                if (maxGenotype.getAllele1HaplotypeId() == 1)
                    isHetHap1.back() = true;
                else if (maxGenotype.getAllele1HaplotypeId() == 2)
                    isHetHap2.back() = true;
            }
        }
    }

    // used to refer isHetHap1, isHetHap2, and numHetVariants
    int bufferActiveRegionIndex = -1;

    // Active region id assigned by the active region detector
    activeRegionId = -1;

    // to record the haplotype id of the first nonref allele of the first variant
    // i.e. to write 0|1 instead of 1|0 for the first variant of the phase block
    uint8_t haplotypeIdOfFirstNonRefAllele(0);

    // to record the first het variant within the phasing set
    pos_t posFirstVariantInPhaseSet(-1);

    // TODO: Comments from Chris
    // Perhaps for a future TDR ticket -- this whole loop looks improvable --
    // it depends on a lot of continuous state between iterations in this loop
    // in ways that could potentially harbor some corner case bugs.
    // I'm not sure how I would suggest to revise,
    // but if we did a cross-team review with other developers,
    // this might be an interesting case to discuss re: refactoring for robustness.
    for (auto& locusPtr : _locusBuffer)
    {
        auto& sampleInfo = locusPtr->getSample(sampleIndex);
        const auto curActiveRegionId(sampleInfo.getActiveRegionId());

        if ((!sampleInfo.isVariant()) || (curActiveRegionId < 0)) continue;

        if (curActiveRegionId != activeRegionId)
        {
            ++bufferActiveRegionIndex;
            posFirstVariantInPhaseSet = -1;

            activeRegionId = curActiveRegionId;
        }

        // no phasing is needed if there's 0 or 1 het variant
        if (numHetVariants[bufferActiveRegionIndex] <= 1)
            continue;

        auto& maxGenotype = sampleInfo.max_gt();

        if ((not maxGenotype.isHet()) or maxGenotype.isConflict()) continue;

        if (posFirstVariantInPhaseSet < 0)
        {
            if (isInstanceOf<GermlineSiteLocusInfo>(*locusPtr))
            {
                // convert to 1-base position
                posFirstVariantInPhaseSet = locusPtr->pos+1;
            }
            else if (isInstanceOf<GermlineIndelLocusInfo>(*locusPtr))
            {
                // for indels no need to add 1 because locusPtr->pos is +1 of the padding base
                posFirstVariantInPhaseSet = locusPtr->pos;
            }
            else
            {
                throw std::bad_cast();
            }
        }

        // phase set id is the 1-based position of the first variant in the set
        sampleInfo.phaseSetId = posFirstVariantInPhaseSet;

        if ((not isHetHap1[bufferActiveRegionIndex]) or (not isHetHap2[bufferActiveRegionIndex]))
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
