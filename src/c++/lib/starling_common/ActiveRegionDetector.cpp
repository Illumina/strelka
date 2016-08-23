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

void ActiveRegionDetector::insertMatch(const align_id_t alignId, const pos_t pos)
{
    setMatch(alignId, pos);
    addAlignIdToPos(alignId, pos);
}

void ActiveRegionDetector::insertSoftClipMatch(const align_id_t alignId, const pos_t pos)
{
    // soft clip doesn't add mismatch count, but the base is used in haplotype generation
    setSoftClipMatch(alignId, pos);
    addAlignIdToPos(alignId, pos);
}

void ActiveRegionDetector::insertSoftClipMismatch(const align_id_t alignId, const pos_t pos, const char baseChar)
{
    // soft clipp doesn't add mismatch count, but the base is used in haplotype generation
    setSoftClipMismatch(alignId, pos, baseChar);
    addAlignIdToPos(alignId, pos);
}

void
ActiveRegionDetector::insertMismatch(const align_id_t alignId, const pos_t pos, const char baseChar)
{
    addCount(pos, MismatchWeight);
    setMismatch(alignId, pos, baseChar);
    addAlignIdToPos(alignId, pos);
}

void
ActiveRegionDetector::insertIndel(const unsigned sampleId, const IndelObservation& indelObservation)
{
    auto pos = indelObservation.key.pos;

    auto alignId = indelObservation.data.id;
    auto indelKey = indelObservation.key;
    if (!indelObservation.data.is_low_map_quality)
    {
        if (indelKey.isPrimitiveInsertionAllele())
        {
            addCount(pos - 1, IndelWeight);
            addCount(pos, IndelWeight);
            setInsert(alignId, pos - 1, indelObservation.key.insert_seq());
            addAlignIdToPos(alignId, pos - 1);
        }
        else if (indelKey.isPrimitiveDeletionAllele())
        {
            unsigned length = indelObservation.key.deletionLength;
            for (unsigned i(0); i<length; ++i)
            {
                addCount(pos + i, IndelWeight);
                setDelete(alignId, pos + i);
                addAlignIdToPos(alignId, pos + i);
            }
            addCount(pos - 1, IndelWeight);
        }
        else
        {
            // ignore BP_LEFT, BP_RIGHT, SWAP
        }
    }
    _indelBuffer.addIndelObservation(sampleId, indelObservation);
}

void
ActiveRegionDetector::updateStartPosition(const pos_t pos)
{
    if (_activeRegions.empty()) return;

    if (_activeRegions.front().getStart() == pos)
    {
        _activeRegions.pop_front();
    }

    _polySites.eraseTo(pos);
}

void
ActiveRegionDetector::updateEndPosition(const pos_t pos, const bool isLastPos)
{
    bool isCurrentPosCandidateVariant = isCandidateVariant(pos);

    // check if we can include this position in the existing acitive region
    bool isSizeFit = (pos - _activeRegionStartPos) < (int)_maxDetectionWindowSize;
    auto distanceFromPrevVariant = pos - _prevVariantPos;
    bool isConsecutiveVariant = (distanceFromPrevVariant == 1); // size may exceed _maxDetectionWindowSize for consecutive variants
    bool isNotFarFromPrevVariant = (pos > 0) && (distanceFromPrevVariant <= MaxDistanceBetweenTwoVariants);

    bool isExtensible = (isSizeFit || isConsecutiveVariant) && isNotFarFromPrevVariant;

    if (isExtensible && !isLastPos)  // if pos is the last position, we cannot extend
    {
        // this position extends the existing active region
        if (isCurrentPosCandidateVariant)
        {
            ++_numVariants;
            _prevVariantPos = pos;
        }
    }
    else
    {
        if (isLastPos && isExtensible && isCurrentPosCandidateVariant)
        {
            ++_numVariants;
            _prevVariantPos = pos;
        }

        // this position doesn't extend the existing active region
        if (_numVariants >= _minNumVariantsPerRegion)
        {
            pos_t origStart = _activeRegionStartPos;
            pos_t origEnd = _prevVariantPos;
            pos_t start, end;
            getExpandedRange(origStart, origEnd, start, end);

            // close existing active region
            std::string refStr = "";
            _ref.get_substring(start, end - start + 1, refStr);

            // TODO: calculate ploidy
            _activeRegions.emplace_back(start, end, refStr, _aligner, _alignIdToAlignInfo);
            ActiveRegion& activeRegion(_activeRegions.back());
            // add haplotype bases
            for (pos_t activeRegionPos(start); activeRegionPos<=end; ++activeRegionPos)
            {
                for (const align_id_t alignId : getPositionToAlignIds(activeRegionPos))
                {
                    std::string haplotypeBase;
                    bool isSoftClipped = setHaplotypeBase(alignId, activeRegionPos, haplotypeBase);
                    if (isSoftClipped)
                        activeRegion.setSoftClipped(alignId);
                    activeRegion.insertHaplotypeBase(alignId, activeRegionPos, haplotypeBase);
                }
            }
            activeRegion.processHaplotypes(_indelBuffer, _polySites);
        }

        if (isCurrentPosCandidateVariant)
        {
            // start new active region
            _activeRegionStartPos = pos;
            _numVariants = 1;
            _prevVariantPos = pos;
        }
        else
        {
            _activeRegionStartPos = 0;
            _numVariants = 0;
        }
    }

    const pos_t minPosWithBase = pos - _maxDetectionWindowSize - _maxDeletionSize;
    clearPos(minPosWithBase);
}

void ActiveRegionDetector::getExpandedRange(const pos_t origStart, const pos_t origEnd, pos_t& newStart, pos_t& newEnd)
{
    // calculate newStart
    newStart = origStart;
    unsigned deltaPos(0);
    for (unsigned repeatUnitLength(1); repeatUnitLength<=MaxRepeatUnitLength; ++repeatUnitLength)
    {
        unsigned repeatSpan = repeatUnitLength;
        for (pos_t pos(origStart-repeatUnitLength); pos >= _ref.get_offset(); --pos)
        {
            char baseChar = _ref.get_base(pos);
            char baseCharToCompare = _ref.get_base(pos+repeatUnitLength);
            if (baseChar != baseCharToCompare)
                break;
            ++repeatSpan;
        }
        unsigned repeatLength = repeatSpan / repeatUnitLength;
        if (repeatLength > 1)
            deltaPos = std::max(deltaPos, repeatSpan);
    }
    deltaPos = std::min(deltaPos, MaxRepeatSpan);
    newStart = origStart - deltaPos;

    // calculate newEnd
    deltaPos = 0;
    for (unsigned repeatUnitLength(1); repeatUnitLength<=MaxRepeatUnitLength; ++repeatUnitLength)
    {
        unsigned repeatSpan = repeatUnitLength;
        for (pos_t pos(origEnd+repeatUnitLength); pos < _ref.end(); ++pos)
        {
            char baseChar = _ref.get_base(pos);
            char baseCharToCompare = _ref.get_base(pos-repeatUnitLength);
            if (baseChar != baseCharToCompare)
                break;
            ++repeatSpan;
        }
        unsigned repeatLength = repeatSpan / repeatUnitLength;
        if (repeatLength > 1)
            deltaPos = std::max(deltaPos, repeatSpan);
    }
    deltaPos = std::min(deltaPos, MaxRepeatSpan);

    newEnd = origEnd + deltaPos;
}

void ActiveRegionDetector::setMatch(const align_id_t id, const pos_t pos)
{
    _variantInfo[id % MaxDepth][pos % MaxBufferSize] = MATCH;
}

void ActiveRegionDetector::setMismatch(const align_id_t id, const pos_t pos, char baseChar)
{
    unsigned idIndex = id % MaxDepth;
    unsigned posIndex = pos % MaxBufferSize;
    _variantInfo[idIndex][posIndex] = MISMATCH;
    _snvBuffer[idIndex][posIndex] = baseChar;
}

void ActiveRegionDetector::setSoftClipMatch(const align_id_t id, const pos_t pos)
{
    _variantInfo[id % MaxDepth][pos % MaxBufferSize] = SOFT_CLIP_MATCH;
}

void ActiveRegionDetector::setSoftClipMismatch(const align_id_t id, const pos_t pos, char baseChar)
{
    unsigned idIndex = id % MaxDepth;
    unsigned posIndex = pos % MaxBufferSize;
    _variantInfo[idIndex][posIndex] = SOFT_CLIP_MISMATCH;
    _snvBuffer[idIndex][posIndex] = baseChar;
}

void ActiveRegionDetector::setDelete(const align_id_t id, const pos_t pos)
{
    _variantInfo[id % MaxDepth][pos % MaxBufferSize] = DELETE;
}

void ActiveRegionDetector::setInsert(const align_id_t id, const pos_t pos, const std::string& insertSeq)
{
    unsigned idIndex = id % MaxDepth;
    unsigned posIndex = pos % MaxBufferSize;
    _variantInfo[idIndex][posIndex] = (_variantInfo[idIndex][posIndex] == MISMATCH ? MISMATCH_INSERT : INSERT);
    _insertSeqBuffer[idIndex][posIndex] = insertSeq;
}

bool ActiveRegionDetector::setHaplotypeBase(const align_id_t id, const pos_t pos, std::string& base) const
{
    unsigned idIndex = id % MaxDepth;
    unsigned posIndex = pos % MaxBufferSize;
    const auto variant = _variantInfo[idIndex][posIndex];
    switch (variant)
    {
    case MATCH:
    case SOFT_CLIP_MATCH:
        base = _ref.get_base(pos);
        break;
    case MISMATCH:
    case SOFT_CLIP_MISMATCH:
        base = std::string(1, _snvBuffer[idIndex][posIndex]);
        break;
    case DELETE:
        base = "";
        break;
    case INSERT:
        base = _ref.get_base(pos) + _insertSeqBuffer[idIndex][posIndex];
        break;
    case MISMATCH_INSERT:
        base = _snvBuffer[idIndex][posIndex] + _insertSeqBuffer[idIndex][posIndex];
    }

    if (variant == SOFT_CLIP_MATCH or variant == SOFT_CLIP_MISMATCH)
        return true;
    return false;
}

bool
ActiveRegionDetector::isCandidateVariant(const pos_t pos) const
{
    auto count = getCount(pos);
    return count >= _minNumVariantsPerPosition && count >= (MinAlternativeAlleleFraction*getDepth(pos));
}

bool ActiveRegionDetector::isPolymorphicSite(const pos_t pos) const
{
    return _polySites.isKeyPresent(pos);
}
