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

void
ActiveRegionDetector::insertMismatch(const align_id_t alignId, const pos_t pos, const char baseChar)
{
    addCount(pos, 1);
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
    _bufferStartPos = pos;

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

    if ((pos - _activeRegionStartPos  >= (int)_maxDetectionWindowSize && pos - _prevVariantPos > 1) || (pos - _prevVariantPos > MaxDistanceBetweenTwoVariants) || isLastPos)
    {
        // this position doesn't extend the existing active region
        if (_numVariants >= _minNumVariantsPerRegion)
        {
            if (isLastPos && isCurrentPosCandidateVariant)
            {
                _prevVariantPos = pos;
            }

            // close existing active region
            std::string refStr = "";
            _ref.get_substring(_activeRegionStartPos, _prevVariantPos - _activeRegionStartPos + 1, refStr);
            _activeRegions.emplace_back(_activeRegionStartPos, _prevVariantPos, refStr, _aligner, _alignIdToAlignInfo);
//            std::cout << _activeRegionStartPos+1 << '\t' << _prevVariantPos+1 << '\t' << refStr << std::endl;
            ActiveRegion& activeRegion(_activeRegions.back());
            // add haplotype bases
            for (pos_t activeRegionPos(_activeRegionStartPos); activeRegionPos<=_prevVariantPos; ++activeRegionPos)
            {
                for (const align_id_t alignId : getPositionToAlignIds(activeRegionPos))
                {
                    std::string haplotypeBase;
                    setHaplotypeBase(alignId, activeRegionPos, haplotypeBase);
                    activeRegion.insertHaplotypeBase(alignId, activeRegionPos, haplotypeBase);
                }
            }
            activeRegion.processHaplotypes(_indelBuffer, _polySites);
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
        // this position extends the existing active region
        if (isCurrentPosCandidateVariant)
        {
            ++_numVariants;
        }
    }

    if (isCurrentPosCandidateVariant)
        _prevVariantPos = pos;

    clearPos(pos - _maxDetectionWindowSize - _maxDeletionSize);
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

void ActiveRegionDetector::setHaplotypeBase(const align_id_t id, const pos_t pos, std::string& base) const
{
    unsigned idIndex = id % MaxDepth;
    unsigned posIndex = pos % MaxBufferSize;
    const auto variant = _variantInfo[idIndex][posIndex];
    switch (variant)
    {
    case MATCH:
        base = _ref.get_base(pos);
        break;
    case MISMATCH:
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
}

//void ActiveRegionDetector::setHaplotypeBaseSnv(const align_id_t id, const pos_t pos, char baseChar)
//{
//    std::string baseStr;
//    switch (baseChar)
//    {
//    case 'A':
//        baseStr = strA;
//        break;
//    case 'C':
//        baseStr = strC;
//        break;
//    case 'G':
//        baseStr = strG;
//        break;
//    case 'T':
//        baseStr = strT;
//        break;
//    }
//    _haplotypeBase[id % MaxDepth][pos % MaxBufferSize] = baseStr;
//}

bool
ActiveRegionDetector::isCandidateVariant(const pos_t pos) const
{
    return getCount(pos) >= _minNumMismatchesPerPosition;
}

bool ActiveRegionDetector::isPolymorphicSite(const pos_t pos) const
{
    return _polySites.isKeyPresent(pos);
}
