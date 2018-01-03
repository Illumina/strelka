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

/// \file
/// \author Sangtae Kim
///

#include "ActiveRegionReadBuffer.hh"

void ActiveRegionReadBuffer::insertMatch(const align_id_t alignId, const pos_t pos)
{
    addVariantCount(pos, 0);
    setMatch(alignId, pos);
    addAlignIdToPos(alignId, pos);
}

void ActiveRegionReadBuffer::insertSoftClipSegment(const align_id_t alignId, const pos_t pos, const std::string& segmentSeq, const bool isBeginEdge)
{
    // For invariant counting
    addVariantCount(pos, IndelWeight);
    if (isBeginEdge)
    {
        addSoftClipCount(pos+1, IndelWeight);
    }
    else
    {
        addSoftClipCount(pos-1, IndelWeight);
    }

    // soft clipp doesn't add mismatch count, but the base is used in haplotype generation
    setSoftClipSegment(alignId, pos, segmentSeq);
    addAlignIdToPos(alignId, pos);
}

void
ActiveRegionReadBuffer::insertMismatch(const align_id_t alignId, const pos_t pos, const char baseChar)
{
    addVariantCount(pos, MismatchWeight);
    setMismatch(alignId, pos, baseChar);
    addAlignIdToPos(alignId, pos);
}

void
ActiveRegionReadBuffer::insertIndel(const IndelObservation& indelObservation)
{
    const auto pos = indelObservation.key.pos;

    const auto alignId = indelObservation.data.id;
    const auto& indelObservationKey = indelObservation.key;
    const auto& indelObservationData = indelObservation.data;
    const auto isExternalCandidate(indelObservationData.is_external_candidate);
    const auto isForcedOutput(indelObservationData.is_forced_output);

    if (!indelObservation.data.is_low_map_quality)
    {
        if (isExternalCandidate and (not isForcedOutput))
        {
            // make sure an active region is created around the external candidate indel
            for (pos_t refPos(pos-1); refPos<indelObservationKey.right_pos(); ++refPos)
            {
                addVariantCount(refPos, MinNumVariantsPerPosition);
            }
            // skip adding this indel to the indel buffer
            return;
        }
        else if (indelObservationKey.isPrimitiveInsertionAllele())
        {
            addVariantCount(pos - 1, IndelWeight);
            addVariantCount(pos, IndelWeight);
            setInsert(alignId, pos - 1, indelObservationKey.insert_seq());
            addAlignIdToPos(alignId, pos - 1);
        }
        else if (indelObservationKey.isPrimitiveDeletionAllele())
        {
            unsigned length = indelObservationKey.deletionLength;
            for (unsigned i(0); i<length; ++i)
            {
                addVariantCount(pos + i, IndelWeight);
                setDelete(alignId, pos + i);
                addAlignIdToPos(alignId, pos + i);
            }
            addVariantCount(pos - 1, IndelWeight);
        }
        else
        {
            // ignore BP_LEFT, BP_RIGHT, SWAP
        }
    }

    _indelBuffer.addIndelObservation(getSampleIndex(alignId), indelObservation);
}

void ActiveRegionReadBuffer::setMatch(const align_id_t id, const pos_t pos)
{
    _variantInfo[id % MaxDepth][pos % MaxBufferSize] = MATCH;
}

void ActiveRegionReadBuffer::setMismatch(const align_id_t id, const pos_t pos, const char baseChar)
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

bool ActiveRegionReadBuffer::getHaplotypeBase(const align_id_t id, const pos_t pos, std::string& base) const
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

    const bool isSoftClipped = (variant == SOFT_CLIP);
    return isSoftClipped;
}

void ActiveRegionReadBuffer::setEndPos(const pos_t endPos)
{
    const bool isPosJumped((_readBufferRange.end_pos()+1) != endPos);

    const auto pos(endPos-1);

    if ((not _readBufferRange.end_pos()) or isPosJumped)
    {
        // initialization
        _readBufferRange.set_begin_pos(pos);
        _refRepeatFinder.initRepeatSpan(pos);
    }

    _refRepeatFinder.updateRepeatSpan(pos + MaxRepeatUnitLength*2u);

    _readBufferRange.set_end_pos(endPos);
}

void ActiveRegionReadBuffer::getReadSegments(
    const known_pos_range2& posRange,
    ActiveRegionReadInfo& readInfo,
    const bool includePartialReads,
    const unsigned minReadSegmentLength) const
{
    std::map<align_id_t, std::string> alignIdToHaplotype;
    std::set<align_id_t> alignIdsReachingEnd;
    std::set<align_id_t> invalidAlignIds;
    std::set<align_id_t> allAlignIds;
    // add haplotype bases
    for (pos_t pos(posRange.begin_pos()); pos<posRange.end_pos(); ++pos)
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

            if (pos == posRange.begin_pos())
                alignIdToHaplotype[alignId] = std::string();

            std::string haplotypeBase;
            const bool isSoftClipped = getHaplotypeBase(alignId, pos, haplotypeBase);
            const bool isContainingN = (haplotypeBase.find('N') != std::string::npos);
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

            if (pos == (posRange.end_pos()-1))    // last position
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
        if (haplotype.length() < minReadSegmentLength) continue;

        readInfo.readSegmentsForHaplotypeGeneration.push_back(entry);
    }

    readInfo.numReadsAlignedToActiveRegion = (unsigned)(allAlignIds.size());
}

bool ActiveRegionReadBuffer::isCandidateVariant(const pos_t pos) const
{
    if (_ref.get_base(pos) == 'N')
        return false;

    const auto count = getVariantCount(pos);
    if ((count >= MinNumVariantsPerPosition and count >= (_minAlternativeAlleleFraction*getDepth(pos)))
        or count >= (MinAlternativeAlleleFractionLowDepth*getDepth(pos)))
        return true;

    return false;
}
