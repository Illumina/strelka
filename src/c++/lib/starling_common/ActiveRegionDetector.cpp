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

void
ActiveRegionDetector::updateStartPosition(const pos_t pos)
{
    for (unsigned sampleId(0); sampleId<_sampleCount; ++sampleId)
    {
        _polySites[sampleId].eraseTo(pos);
    }
    _readBuffer.clearPos(pos);
}

void
ActiveRegionDetector::processActiveRegion()
{
    if (_activeRegionPtr)
    {
        _activeRegionPtr->processHaplotypes(_indelBuffer, _polySites);
        _activeRegionPtr.reset();
    }
}

void
ActiveRegionDetector::updateEndPosition(const pos_t pos)
{
    _readBuffer.setEndPos(pos+1);

    const pos_t posToClear = pos - _maxDetectionWindowSize - _maxIndelSize;

    if (pos < _lastActiveRegionEnd)
    {
        // this position was already covered by the previous active region
        _activeRegionStartPos = 0;
        _numVariants = 0;
        _readBuffer.clearPos(posToClear);
        return;
    }

    pos_t posToProcess(pos-1);

    bool isCurrentPosCandidateVariant = isCandidateVariant(posToProcess);

    // check if we can include this position in the existing acitive region
    bool isSizeFit = (posToProcess - _activeRegionStartPos) < (int)_maxDetectionWindowSize;
    auto distanceFromPrevVariant = posToProcess - _prevVariantPos;
    bool isConsecutiveVariant = (distanceFromPrevVariant == 1); // size may exceed _maxDetectionWindowSize for consecutive variants
    bool isNotFarFromPrevVariant = (posToProcess > 0) && (distanceFromPrevVariant <= MaxDistanceBetweenTwoVariants);

    bool isExtensible = (isSizeFit || isConsecutiveVariant) && isNotFarFromPrevVariant;

    if (isExtensible)  // if pos is the last position, we cannot extend
    {
        // this position extends the existing active region
        if (isCurrentPosCandidateVariant)
        {
            ++_numVariants;
            _prevVariantPos = posToProcess;
        }
    }
    else
    {
        // this position doesn't extend the existing active region
        if (_numVariants >= _minNumVariantsPerRegion)
        {
            pos_t origBeginPos = _activeRegionStartPos;
            pos_t origEndPos = _prevVariantPos + 1;
            pos_range newActiveRegion;

            // expand active region to include repeats
            getExpandedRange(pos_range(origBeginPos, origEndPos), newActiveRegion);
            _lastActiveRegionEnd = newActiveRegion.end_pos;

            if (_activeRegionPtr and (_activeRegionPtr->getEndPosition() >= newActiveRegion.begin_pos))
            {
                // merge with the previous active region
                _activeRegionPtr->extendEndPosition(newActiveRegion.end_pos);
            }
            else
            {
                _activeRegionPtr = std::unique_ptr<ActiveRegion>(
                        new ActiveRegion(newActiveRegion, _ref, _maxIndelSize, _sampleCount,
                                         _aligner, _alignerForAssembly, _readBuffer));
            }
        }

        if (isCurrentPosCandidateVariant)
        {
            // start new active region
            _activeRegionStartPos = posToProcess;
            _numVariants = 1;
            _prevVariantPos = posToProcess;
        }
        else
        {
            _activeRegionStartPos = 0;
            _numVariants = 0;
        }
    }

    if (_activeRegionPtr and (_activeRegionPtr->getEndPosition() <= pos))
    {
        processActiveRegion();
    }
}

void ActiveRegionDetector::clear()
{
    if (_numVariants >= _minNumVariantsPerRegion)
    {
        pos_t origBeginPos = _activeRegionStartPos;
        pos_t origEndPos = _prevVariantPos + 1;
        pos_range newActiveRegion;

        // expand active region to include repeats
        getExpandedRange(pos_range(origBeginPos, origEndPos), newActiveRegion);
        _lastActiveRegionEnd = newActiveRegion.end_pos;

        if (_activeRegionPtr and (_activeRegionPtr->getEndPosition() >= newActiveRegion.begin_pos))
        {
            // merge with the previous active region
            _activeRegionPtr->extendEndPosition(newActiveRegion.end_pos);
        }
        else
        {
            _activeRegionPtr = std::unique_ptr<ActiveRegion>(
                    new ActiveRegion(newActiveRegion, _ref, _maxIndelSize, _sampleCount,
                                     _aligner, _alignerForAssembly, _readBuffer));
        }
    }

    processActiveRegion();

    _numVariants = 0;
    _activeRegionStartPos = 0;
    _prevVariantPos = 0;
    _lastActiveRegionEnd = 0;
}

void ActiveRegionDetector::getExpandedRange(const pos_range& origActiveRegion, pos_range& newActiveRegion)
{
    // if origStart is within repeat region, move the start position outside of the repeat region
    // e.g. In TACGACGAC|GAC if origStart points C before |, the start position is moved to T
    // note that bases after origStart are ignored. For example, in TACGAC|GAC, the start position doesn't move.
    pos_t origStart = origActiveRegion.begin_pos;
    const pos_t minStartLowerBound(std::max(0,_lastActiveRegionEnd-1));
    const pos_t minStart(std::max(origStart - MaxRepeatSpan, minStartLowerBound));
    pos_t newBeginPos;
    for (newBeginPos = origStart; newBeginPos > minStart; --newBeginPos)
    {
        bool isPrevPosLowDepth = false;
        for (unsigned sampleId(0); sampleId<_sampleCount; ++sampleId)
        {
            if (_readBuffer.getDepth(sampleId, newBeginPos-1) < MinDepth)
            {
                isPrevPosLowDepth = true;
                break;
            }
        }
        if (isPrevPosLowDepth) break;

        bool isAnchor = _readBuffer.isAnchor(newBeginPos) and (not isCandidateVariant(newBeginPos));
        if (isAnchor) break;
    }
    newActiveRegion.set_begin_pos(newBeginPos);

    // calculate newEnd
    pos_t origEnd = origActiveRegion.end_pos;
    pos_t maxEnd(std::min(origEnd + MaxRepeatSpan, _ref.end()));
    pos_t newEndPos;
    for (newEndPos = origEnd; newEndPos < maxEnd; ++newEndPos)
    {
        bool isPrevPosLowDepth = false;
        for (unsigned sampleId(0); sampleId<_sampleCount; ++sampleId)
        {
            if (_readBuffer.getDepth(sampleId, newEndPos-1) < MinDepth)
            {
                isPrevPosLowDepth = true;
                break;
            }
        }
        if (isPrevPosLowDepth) break;
        bool isAnchor = _readBuffer.isAnchor(newEndPos) and (not isCandidateVariant(newEndPos));
        if (isAnchor) break;
    }
    newActiveRegion.set_end_pos(newEndPos);
}

bool
ActiveRegionDetector::isCandidateVariant(const pos_t pos) const
{
    for (unsigned sampleId(0); sampleId<_sampleCount; ++sampleId)
    {
        if (_ref.get_base(pos) == 'N')
            return false;
        auto count = _readBuffer.getVariantCount(sampleId, pos);
        if ((count >= _minNumVariantsPerPositionPerSample and count >= (MinAlternativeAlleleFraction*_readBuffer.getDepth(sampleId, pos)))
            or count >= (0.4*_readBuffer.getDepth(sampleId, pos)))
            return true;
    }
    return false;
}

bool ActiveRegionDetector::isPolymorphicSite(const unsigned sampleId, const pos_t pos) const
{
    return _polySites[sampleId].isKeyPresent(pos);
}
