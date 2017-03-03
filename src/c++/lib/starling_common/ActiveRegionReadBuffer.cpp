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

#include "ActiveRegionReadBuffer.hh"

void ActiveRegionReadBuffer::insertMatch(const align_id_t alignId, const pos_t pos)
{
    addVariantCount(getSampleId(alignId), pos, 0);
    setMatch(alignId, pos);
    addAlignIdToPos(alignId, pos);
}

void ActiveRegionReadBuffer::insertSoftClipSegment(const align_id_t alignId, const pos_t pos, const std::string& segmentSeq, bool isBeginEdge)
{
    // For invariant counting
    addVariantCount(getSampleId(alignId), pos, IndelWeight);
    unsigned sampleId(getSampleId(alignId));
    if (isBeginEdge)
    {
        addSoftClipCount(sampleId, pos+1, IndelWeight);
    }
    else
    {
        addSoftClipCount(sampleId, pos-1, IndelWeight);
    }

    // soft clipp doesn't add mismatch count, but the base is used in haplotype generation
    setSoftClipSegment(alignId, pos, segmentSeq);
    addAlignIdToPos(alignId, pos);
}

void
ActiveRegionReadBuffer::insertMismatch(const align_id_t alignId, const pos_t pos, const char baseChar)
{
    addVariantCount(getSampleId(alignId), pos, MismatchWeight);
    setMismatch(alignId, pos, baseChar);
    addAlignIdToPos(alignId, pos);
}

void
ActiveRegionReadBuffer::insertIndel(const unsigned sampleId, const IndelObservation& indelObservation)
{
    auto pos = indelObservation.key.pos;

    auto alignId = indelObservation.data.id;
    auto indelKey = indelObservation.key;

    if (!indelObservation.data.is_low_map_quality)
    {
        if (indelKey.isPrimitiveInsertionAllele())
        {
            addVariantCount(sampleId, pos - 1, IndelWeight);
            addVariantCount(sampleId, pos, IndelWeight);
            setInsert(alignId, pos - 1, indelObservation.key.insert_seq());
            addAlignIdToPos(alignId, pos - 1);
        }
        else if (indelKey.isPrimitiveDeletionAllele())
        {
            unsigned length = indelObservation.key.deletionLength;
            for (unsigned i(0); i<length; ++i)
            {
                addVariantCount(sampleId, pos + i, IndelWeight);
                setDelete(alignId, pos + i);
                addAlignIdToPos(alignId, pos + i);
            }
            addVariantCount(sampleId, pos - 1, IndelWeight);
        }
        else
        {
            // ignore BP_LEFT, BP_RIGHT, SWAP
        }
    }
    _indelBuffer.addIndelObservation(sampleId, indelObservation);
}

void ActiveRegionReadBuffer::setMatch(const align_id_t id, const pos_t pos)
{
    _variantInfo[id % MaxDepth][pos % MaxBufferSize] = MATCH;
}

void ActiveRegionReadBuffer::setMismatch(const align_id_t id, const pos_t pos, char baseChar)
{
    unsigned idIndex = id % MaxDepth;
    unsigned posIndex = pos % MaxBufferSize;
    _variantInfo[idIndex][posIndex] = MISMATCH;
    _snvBuffer[idIndex][posIndex] = baseChar;
}

void ActiveRegionReadBuffer::setSoftClipSegment(const align_id_t id, const pos_t pos, const std::string& segmentSeq)
{
    unsigned idIndex = id % MaxDepth;
    unsigned posIndex = pos % MaxBufferSize;
    _variantInfo[idIndex][posIndex] = SOFT_CLIP;
    _insertSeqBuffer[idIndex][posIndex] = segmentSeq;
}

void ActiveRegionReadBuffer::setDelete(const align_id_t id, const pos_t pos)
{
    _variantInfo[id % MaxDepth][pos % MaxBufferSize] = DELETE;
}

void ActiveRegionReadBuffer::setInsert(const align_id_t id, const pos_t pos, const std::string& insertSeq)
{
    unsigned idIndex = id % MaxDepth;
    unsigned posIndex = pos % MaxBufferSize;
    _variantInfo[idIndex][posIndex] = (_variantInfo[idIndex][posIndex] == MISMATCH ? MISMATCH_INSERT : INSERT);
    _insertSeqBuffer[idIndex][posIndex] = insertSeq;
}

bool ActiveRegionReadBuffer::getHaplotypeBase(const align_id_t id, const pos_t pos, std::string &base) const
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
    case SOFT_CLIP:
        base = _insertSeqBuffer[idIndex][posIndex];
        break;
    case MISMATCH_INSERT:
        base = _snvBuffer[idIndex][posIndex] + _insertSeqBuffer[idIndex][posIndex];
    }

    bool isSoftClipped = (variant == SOFT_CLIP);
    return isSoftClipped;
}

void ActiveRegionReadBuffer::setEndPos(pos_t endPos)
{
    auto pos(endPos-1);

    if (not _readBufferRange.is_begin_pos)
    {
        // initialization
        _readBufferRange.set_begin_pos(pos);
        _refRepeatFinder.initRepeatSpan(pos);
    }

    _refRepeatFinder.updateRepeatSpan(pos + MaxRepeatUnitLength*2u);

    _readBufferRange.set_end_pos(endPos);
}

void ActiveRegionReadBuffer::getReadSegments(pos_range posRange, ReadInfo& readInfo, bool includePartialReads) const
{
    std::map<align_id_t, std::string> alignIdToHaplotype;
    std::set<align_id_t> alignIdsReachingEnd;
    std::set<align_id_t> invalidAlignIds;
    std::set<align_id_t> allAlignIds;
    // add haplotype bases
    for (pos_t pos(posRange.begin_pos); pos<posRange.end_pos; ++pos)
    {
        if (not _readBufferRange.is_pos_intersect(pos))
            continue;

        for (const align_id_t alignId : getPositionToAlignIds(pos))
        {
            allAlignIds.insert(alignId);
            if (not includePartialReads)
            {
                if (invalidAlignIds.count(alignId)) continue;
            }

            if (pos == posRange.begin_pos)
                alignIdToHaplotype[alignId] = std::string();

            std::string haplotypeBase;
            bool isSoftClipped = getHaplotypeBase(alignId, pos, haplotypeBase);
            bool isContainingN = (haplotypeBase.find('N') != std::string::npos);
            if (isContainingN or (!includePartialReads and isSoftClipped))
                invalidAlignIds.insert(alignId);

            if (!alignIdToHaplotype.count(alignId))
            {
                if (includePartialReads)
                    alignIdToHaplotype[alignId] = std::string();
                else
                    continue;   // this read doesn't cover beginPos, ignore it
            }
            if (includePartialReads or (not isSoftClipped))
                alignIdToHaplotype[alignId] += haplotypeBase;

            if (pos == (posRange.end_pos-1))    // last position
                alignIdsReachingEnd.insert(alignId);
        }
    }

    for (const auto& entry : alignIdToHaplotype)
    {
        align_id_t alignId = entry.first;

        if (!includePartialReads)
        {
            if (invalidAlignIds.count(alignId) or !alignIdsReachingEnd.count(alignId))
                continue;
        }

        const auto& haplotype(entry.second);
        if (haplotype.empty()) continue;

        readInfo.readSegments.push_back(entry);
    }

    readInfo.numReads = (unsigned)(allAlignIds.size());
}

bool ActiveRegionReadBuffer::isCandidateVariant(const pos_t pos) const
{
    for (unsigned sampleId(0); sampleId<_sampleCount; ++sampleId)
    {
        if (_ref.get_base(pos) == 'N')
            return false;
        auto count = getVariantCount(sampleId, pos);
        if ((count >= MinNumVariantsPerPosition and count >= (MinAlternativeAlleleFraction*getDepth(sampleId, pos)))
            or count >= (MinAlternativeAlleleFractionLowDepth*getDepth(sampleId, pos)))
            return true;
    }
    return false;
}