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

///
/// \author Sangtae Kim
///

#include "ActiveRegionDetector.hh"
#include <iostream>

void ActiveRegionDetector::insertMatch(const align_id_t alignId, const pos_t pos, const char baseCode)
{
    setHaplotypeBaseSnv(alignId, pos, baseCode);
    addAlignIdToPos(alignId, pos);
}

void
ActiveRegionDetector::insertMismatch(const align_id_t alignId, const pos_t pos, const char baseChar)
{
    addCount(pos, 1);
    setHaplotypeBaseSnv(alignId, pos, baseChar);
    addAlignIdToPos(alignId, pos);
}

void
ActiveRegionDetector::insertIndel(const IndelObservation& indelObservation)
{
    auto pos = indelObservation.key.pos;
    const int indelCount = 4;

    auto alignId = indelObservation.data.id;
    auto indelKey = indelObservation.key;
    if (indelKey.type == INDEL::INSERT)
    {
        addCount(pos - 1, indelCount);
        addCount(pos, indelCount);
        concatenateHaplotypeBase(alignId, pos - 1, indelObservation.data.insert_seq);
        addAlignIdToPos(alignId, pos - 1);
    }
    else if (indelKey.type == INDEL::DELETE)
    {
        unsigned length = indelObservation.key.length;
        for (unsigned i(0); i<length; ++i)
        {
            addCount(pos + i, indelCount);
            setHaplotypeBase(alignId, pos + i, "");
            addAlignIdToPos(alignId, pos + i);
        }
        addCount(pos - 1, indelCount);
    }
    else
    {
        // ignore BP_LEFT, BP_RIGHT, SWAP
    }
}

//static int activeRegionIndex = 0;

void
ActiveRegionDetector::updateStartPosition(const pos_t pos)
{
    for (int i(_bufferStartPos); i<pos; ++i) resetCounter(i);

    _bufferStartPos = pos;

    if (_activeRegions.empty()) return;

    ActiveRegion front =  _activeRegions.front();
    if (front.getStart() == pos)
    {
        _activeRegions.pop_front();
    }
}

void
ActiveRegionDetector::updateEndPosition(const pos_t pos)
{
    bool isCurrentPosCandidateVariant = isCandidateVariant(pos);
    if (pos - _activeRegionStartPos  >= (int)_maxDetectionWindowSize && pos - _prevVariantPos > 1)
    {
        // this position doesn't extend the existing active region
        if (_numVariants >= _minNumVariantsPerRegion)
        {
            // close existing active region
            ActiveRegion activeRegion = ActiveRegion(_activeRegionStartPos, _prevVariantPos, _numVariants);
            _activeRegions.push_back(activeRegion);
            // add haplotype bases
            std::cout << '>' << (_activeRegionStartPos+1) << '\t' << (_prevVariantPos+1) << '\t' << _numVariants << std::endl;
            for (pos_t activeRegionPos(_activeRegionStartPos); activeRegionPos<=_prevVariantPos; ++activeRegionPos)
            {
                for (align_id_t alignId : getPositionToAlignIds(activeRegionPos))
                {
                    activeRegion.insertHaplotypeBase(alignId, activeRegionPos,
                                                     getHaplotypeBase(alignId, activeRegionPos));
                }
            }
            activeRegion.printHaplotypeSequences();
        }
        if (!isCurrentPosCandidateVariant)
        {
            _activeRegionStartPos = -1;
            _numVariants = 0;
        }
        else
        {
            // start new active region
            _activeRegionStartPos = pos;
            _numVariants = 1;
        }
    }
    else
    {
        // this position doesn't extend the existing active region
        if (isCurrentPosCandidateVariant)
        {
            ++_numVariants;
        }
    }

    if (isCurrentPosCandidateVariant)
        _prevVariantPos = pos;

    clearPos(pos - _maxDetectionWindowSize - 50);

    // max_deletion_size = 50
    // TODO: remaining
}

void ActiveRegionDetector::setHaplotypeBaseSnv(const align_id_t id, const pos_t pos, char baseChar)
{
    std::string baseStr;
    switch (baseChar)
    {
    case 'A':
        baseStr = strA;
        break;
    case 'C':
        baseStr = strC;
        break;
    case 'G':
        baseStr = strG;
        break;
    case 'T':
        baseStr = strT;
        break;
    }
    _haplotypeBase[id % MaxDepth][pos % MaxBufferSize] = baseStr;
}

bool
ActiveRegionDetector::isCandidateVariant(const pos_t pos) const
{
    return getCount(pos) >= _minNumMismatchesPerPosition;
}

std::list<ActiveRegion>&
ActiveRegionDetector::getActiveRegions()
{
    return _activeRegions;
}
