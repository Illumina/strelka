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

/// \file
/// \author Sangtae Kim
///

#include "ActiveRegionDetector.hh"


void
ActiveRegionDetector::clearReadBuffer(const pos_t pos)
{
    _readBuffer.clearPos(pos);
}

void
ActiveRegionDetector::clearUpToPos(const pos_t pos)
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
ActiveRegionDetector::closeExistingActiveRegion()
{
    // close the existing active region
    assert (_activeRegionStartPos < _anchorPosFollowingPrevVariant);

    pos_range activeRegionRange(_activeRegionStartPos, _anchorPosFollowingPrevVariant + 1);
    _activeRegions.emplace_back(activeRegionRange, _ref, _maxIndelSize, _sampleIndex,
                                _aligner, _readBuffer, _indelBuffer, _candidateSnvBuffer);
    setPosToActiveRegionIdMap(activeRegionRange);

    // we have no existing active region at this point
    _numVariants = 0;
    _activeRegionStartPos = 0;
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
    const bool isDepthZero = _readBuffer.isDepthZero(posToProcess);
    // depth 0 position can be a candidate variant pos
    // but cannot open a new active region
    if (isDepthZero and (_numVariants == 0u))
    {
        isCurrentPosCandidateVariant = false;
    }

    const bool isAnchor = _readBuffer.isAnchor(posToProcess) and (not isCurrentPosCandidateVariant);

    if (!isCurrentPosCandidateVariant and !isAnchor) return;

//    std::cout << (posToProcess+1) << "\t" << isCurrentPosCandidateVariant << '\t' << isAnchor << std::endl;

    unsigned distanceFromPrevVariant = (unsigned) (posToProcess - _prevVariantPos);
    if (distanceFromPrevVariant > MaxDistanceBetweenTwoVariants and _anchorPosFollowingPrevVariant >= 0)
    {
        if (_numVariants >= MinNumVariantsPerRegion )
        {
            closeExistingActiveRegion();
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

        if (_anchorPosFollowingPrevVariant < 0)
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
        _anchorPosFollowingPrevVariant = -1;
    }
}

void ActiveRegionDetector::clear()
{
    if (_isBeginning) return;

    if (_numVariants >= MinNumVariantsPerRegion)
    {
        bool isStartPosUndetermined = (_activeRegionStartPos < 0);
        bool isEndPosUndetermined = (_anchorPosFollowingPrevVariant < 0);

        // if both the start and end pos are undetermined
        // and read buffer size is zero, do not close the active region
        if (!isStartPosUndetermined || !isEndPosUndetermined
            || (_readBuffer.getBeginPos() < _readBuffer.getEndPos()))
        {
            if (isStartPosUndetermined)
                _activeRegionStartPos = _readBuffer.getBeginPos();
            if (isEndPosUndetermined)
                _anchorPosFollowingPrevVariant = (_readBuffer.getEndPos()-1);
            closeExistingActiveRegion();
        }
    }

    processActiveRegion();

    _activeRegionStartPos = -1;
    _anchorPosFollowingPrevVariant = -1;
    _prevVariantPos = -1;
    _numVariants = 0;
}
