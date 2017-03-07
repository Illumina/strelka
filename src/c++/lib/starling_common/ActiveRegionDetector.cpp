// -*- mode: c++; indent-tabs-mode: nil; -*-
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

///
/// \author Sangtae Kim
///

#include "ActiveRegionDetector.hh"

void
ActiveRegionDetector::clearReadBuffer(const pos_t pos)
{
    _readBuffer.clearPos(pos);
}

void
ActiveRegionDetector::clearPolySites(const pos_t pos)
{
    for (unsigned sampleId(0); sampleId<_sampleCount; ++sampleId)
    {
        _polySites[sampleId].eraseTo(pos);
    }
}

void
ActiveRegionDetector::clearPosToActiveRegionMap(const pos_t pos)
{
    _posToActiveRegionIdMap.eraseTo(pos);
}


void
ActiveRegionDetector::setPosToActiveRegionIdMap(pos_range activeRegionRange)
{
    if (activeRegionRange.size() > ActiveRegion::MaxRefSpanToBypassAssembly)
        return;

    ActiveRegionId activeRegionId(activeRegionRange.begin_pos);
    for (pos_t pos(activeRegionRange.begin_pos); pos<activeRegionRange.end_pos; ++pos)
    {
        _posToActiveRegionIdMap.getRef(pos) = activeRegionId;
    }
}



void
ActiveRegionDetector::processActiveRegion()
{
    if (not _activeRegions.empty())
    {
        _activeRegions.front().processHaplotypes();
        _activeRegions.pop_front();
    }
}


void
ActiveRegionDetector::updateEndPosition(const pos_t pos)
{
    _readBuffer.setEndPos(pos+1);

    if (_isBeginning)
    {
        _activeRegionStartPos = pos;
        _anchorPosFollowingPrevVariant = pos;
        _prevAnchorPos = pos;
        _isBeginning = false;
    }

    // process and pop active regions
    if ((not _activeRegions.empty()) and (_activeRegions.front().getEndPosition() <= (pos_t)(pos-ActiveRegionReadBuffer::MaxAssemblyPadding)))
    {
        processActiveRegion();
    }

    pos_t posToProcess(pos-1);
    if (posToProcess < 0) return;

    bool isCurrentPosCandidateVariant = _readBuffer.isCandidateVariant(posToProcess);
    bool isAnchor = _readBuffer.isAnchor(posToProcess) and (not isCurrentPosCandidateVariant);

    if (!isCurrentPosCandidateVariant and !isAnchor) return;

//    std::cout << (posToProcess+1) << "\t" << isCurrentPosCandidateVariant << '\t' << isAnchor << std::endl;

    unsigned distanceFromPrevVariant = (unsigned) (posToProcess - _prevVariantPos);
    if (distanceFromPrevVariant > MaxDistanceBetweenTwoVariants and _anchorPosFollowingPrevVariant > 0)
    {
        if (_numVariants >= MinNumVariantsPerRegion )
        {
            // close the existing active region
            pos_range activeRegionRange(_activeRegionStartPos, _anchorPosFollowingPrevVariant + 1);
            _activeRegions.emplace_back(activeRegionRange, _ref, _maxIndelSize, _sampleCount,
                                        _aligner, _readBuffer, _indelBuffer, _polySites);

            setPosToActiveRegionIdMap(activeRegionRange);

            // we have no existing acive region at this point
            _numVariants = 0;
            _activeRegionStartPos = 0;
        }
        else
        {
            // no active region is detected
            _numVariants = 0;
        }
    }

    if (isAnchor)
    {
        if (_numVariants == 0)
        {
            // this position can be a start of an active region
            _activeRegionStartPos = posToProcess;
        }

        if (not _anchorPosFollowingPrevVariant)
        {
            _anchorPosFollowingPrevVariant = posToProcess;
        }
        _prevAnchorPos = posToProcess;
    }

    if (isCurrentPosCandidateVariant)
    {
        if (not _activeRegionStartPos)
            _activeRegionStartPos = _prevAnchorPos;
        // extend the existing active region
        ++_numVariants;
        _prevVariantPos = posToProcess;
        _anchorPosFollowingPrevVariant = 0;
    }
}

void ActiveRegionDetector::clear()
{
    if (_isBeginning) return;

    if (_numVariants >= MinNumVariantsPerRegion)
    {
        // close the existing active region
        if (not _anchorPosFollowingPrevVariant)
            _anchorPosFollowingPrevVariant = _readBuffer.getEndPos();
        pos_range activeRegionRange(_activeRegionStartPos, _anchorPosFollowingPrevVariant + 1);
        _activeRegions.emplace_back(activeRegionRange, _ref, _maxIndelSize, _sampleCount,
                                    _aligner, _readBuffer, _indelBuffer, _polySites);
        setPosToActiveRegionIdMap(activeRegionRange);
    }

    processActiveRegion();

    _activeRegionStartPos = 0;
    _anchorPosFollowingPrevVariant = 1;
    _prevVariantPos = 0;
    _numVariants = 0;
}

bool ActiveRegionDetector::isPolymorphicSite(const unsigned sampleId, const pos_t pos) const
{
    return _polySites[sampleId].isKeyPresent(pos);
}

uint8_t ActiveRegionDetector::getHaplotypeId(const unsigned sampleId, const pos_t pos, const BASE_ID::index_t baseIndex) const
{
    if (not isPolymorphicSite(sampleId, pos)) return 0; // reference (complex allele index 0)
    auto value(_polySites[sampleId].getConstRef(pos));
    return ActiveRegion::getHaplotypeId(value, baseIndex);
}