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

#pragma once

#include <blt_util/reference_contig_segment.hh>
#include <vector>
#include <iostream>
#include "starling_types.hh"
#include "indel.hh"
#include "IndelBuffer.hh"
#include "ReferenceRepeatFinder.hh"

/// AlignInfo object to store sample id and indel align type
struct AlignInfo
{
    unsigned sampleId;
    INDEL_ALIGN_TYPE::index_t indelAlignType;
};

struct ReadInfo
{
    std::vector<std::pair<align_id_t, std::string>> readSegments;
    unsigned numReads;
};

class ActiveRegionReadBuffer
{
public:

    // maximum buffer size in bases (must be larger than the maximum read size + max indel size
    static const unsigned MaxBufferSize = 1000u;

    // maximum read depth
    // TODO: dynamically calculate maximum depth
    static const unsigned MaxDepth = 1000u;

    static const unsigned MinNumVariantsPerPosition = 9u;

    static const pos_t MaxAssemblyPadding = (pos_t)9u;

    // maximum repeat unit to consider
    static const unsigned MaxRepeatUnitLength = 50u;

    static const unsigned MinRepeatSpan = (pos_t)3u;

    // variant count to add for a single mismatch or indel
    static const int MismatchWeight = 1;
    static const int IndelWeight = 4;

    // minimum alternative allele fraction to call a position as a candidate variant
    const float MinAlternativeAlleleFraction = 0.2;

    // if the fraction is larger than MinAlternativeAlleleFractionLowDepth
    // the position becomes candidate even if the number is lower than MinNumVariantsPerPosition
    const float MinAlternativeAlleleFractionLowDepth = 0.35;

    /// Read buffer to be used in active regions
    /// \param ref reference
    /// \param sampleCount sample count
    /// \param indelBuffer indel buffer
    ActiveRegionReadBuffer(
            const reference_contig_segment& ref,
            unsigned sampleCount,
            IndelBuffer& indelBuffer)
            :
            _ref(ref),
            _refRepeatFinder(ref, MaxRepeatUnitLength, MaxBufferSize, MinRepeatSpan),
            _sampleCount(sampleCount),
            _indelBuffer(indelBuffer),
            _variantCounter(sampleCount, std::vector<unsigned>(MaxBufferSize)),
            _depth(sampleCount, std::vector<unsigned>(MaxBufferSize)),
            _positionToAlignIds(MaxBufferSize),
            _alignIdToAlignInfo(MaxDepth),
            _variantInfo(MaxDepth, std::vector<VariantType>(MaxBufferSize, VariantType())),
            _insertSeqBuffer(MaxDepth, std::vector<std::string>(MaxBufferSize, std::string()))
    {}

    /// insert match at position pos
    /// \param alignId align id
    /// \param pos reference position
    void insertMatch(const align_id_t alignId, const pos_t pos);

    /// insert mismatch at position pos
    /// \param alignId align id
    /// \param pos reference position
    /// \param baseChar read base char
    void insertMismatch(const align_id_t alignId, const pos_t pos, const char baseChar);

    /// insert soft-clipped segment
    /// \param alignId align id
    /// \param pos reference position
    /// \param baseChar soft-clipped segment sequence
    void insertSoftClipSegment(const align_id_t alignId, const pos_t pos, const std::string& segmentSeq, bool isBeginEdge);

    /// insert indel
    /// \param sampleId sample id
    /// \param indelObservation indel observation object
    void insertIndel(const unsigned sampleId, const IndelObservation& indelObservation);

    /// checks if pos is an anchor position
    /// \param pos reference position
    /// \return true if pos is an anchor position, false otherwise
    bool isAnchor(pos_t pos) const
    {
        return _refRepeatFinder.isAnchor(pos);
    }

    /// Set end position of the buffer
    /// \param endPos end position (exclusive)
    void setEndPos(pos_t endPos);

    /// Gets the beginning position
    /// \return begin position
    pos_t getBeginPos() const { return _readBufferRange.begin_pos; }

    /// Gets the end position
    /// \return end position
    pos_t getEndPos() const { return _readBufferRange.end_pos; }

    /// Gets sample id and indel align type
    /// \param alignId align id
    /// \return sample id and indel align type
    const AlignInfo& getAlignInfo(align_id_t alignId) const
    {
        return _alignIdToAlignInfo[alignId % MaxDepth];
    }

    /// Gets read segments for the range
    /// \param posRange position range
    /// \param readInfo read info object to store read segments
    /// \param includePartialReads if true, only reads fully covering the region will be retrieved
    void getReadSegments(pos_range posRange, ReadInfo &readInfo, bool includePartialReads) const;

    /// cache sampleId and indelAlignType corresponding to alignId
    /// \param alignId align id
    /// \param sampleId sample id
    /// \param indelAlignType indel align type
    void setAlignInfo(const align_id_t alignId, unsigned sampleId, INDEL_ALIGN_TYPE::index_t indelAlignType)
    {
        AlignInfo& alignInfo = _alignIdToAlignInfo[alignId % MaxDepth];
        alignInfo.sampleId = sampleId;
        alignInfo.indelAlignType = indelAlignType;
    }

    /// Gets sample id
    /// \param alignId align id
    /// \return sample id
    unsigned getSampleId(const align_id_t alignId) const
    {
        return _alignIdToAlignInfo[alignId % MaxDepth].sampleId;
    }

    /// Clear buffer
    /// \param pos position
    void clearPos(pos_t pos)
    {
        _positionToAlignIds[pos % MaxBufferSize].clear();
        resetCounter(pos);
        _readBufferRange.set_begin_pos(pos+1);
    }

    /// Checks if pos is candidate variant position
    /// \param pos position
    /// \return true if pos is candidate variant, false otherwise
    bool isCandidateVariant(const pos_t pos) const;

private:
    enum VariantType
    {
        MATCH,
        MISMATCH,
        SOFT_CLIP,
        DELETE,
        INSERT,
        MISMATCH_INSERT
    };

    const reference_contig_segment& _ref;
    ReferenceRepeatFinder _refRepeatFinder;

    const unsigned _sampleCount;
    IndelBuffer& _indelBuffer;

    pos_range _readBufferRange;

    std::vector<std::vector<unsigned>> _variantCounter;
    std::vector<std::vector<unsigned>> _depth;

    // for haplotypes
    std::vector<std::vector<align_id_t>> _positionToAlignIds;

    // to store align information
    std::vector<AlignInfo> _alignIdToAlignInfo;

    std::vector<std::vector<VariantType>> _variantInfo;
    std::vector<std::vector<std::string>> _insertSeqBuffer;
    char _snvBuffer[MaxDepth][MaxBufferSize];

    void setMatch(const align_id_t id, const pos_t pos);
    void setMismatch(const align_id_t id, const pos_t pos, char baseChar);
    void setDelete(const align_id_t id, const pos_t pos);
    void setInsert(const align_id_t id, const pos_t pos, const std::string& insertSeq);
    void setSoftClipSegment(const align_id_t id, const pos_t pos, const std::string& segmentSeq);
    bool setHaplotypeBase(const align_id_t id, const pos_t pos, std::string& base) const;

    void resetCounter(const pos_t pos)
    {
        int index = pos % MaxBufferSize;

        for (unsigned sampleId=0; sampleId<_sampleCount; ++sampleId)
        {
            _variantCounter[sampleId][index] = 0;
            _depth[sampleId][index] = 0;
        }
    }

    void addVariantCount(const unsigned sampleId, const pos_t pos, unsigned count)
    {
        int index = pos % MaxBufferSize;
        _variantCounter[sampleId][index] += count;
        ++_depth[sampleId][index];
    }

    void addSoftClipCount(const unsigned sampleId, const pos_t pos, unsigned count)
    {
        int index = pos % MaxBufferSize;
        _variantCounter[sampleId][index] += count;
    }

    unsigned getVariantCount(const unsigned sampleId, const pos_t pos) const
    {
        return _variantCounter[sampleId][pos % MaxBufferSize];
    }

    void addAlignIdToPos(const align_id_t alignId, const pos_t pos)
    {
        int index = pos % MaxBufferSize;
        if (_positionToAlignIds[index].empty() || _positionToAlignIds[index].back() != alignId)
            _positionToAlignIds[index].push_back(alignId);
    }

    unsigned getDepth(const unsigned sampleId, const pos_t pos) const
    {
        int index = pos % MaxBufferSize;
        return _depth[sampleId][index];
    }

    const std::vector<align_id_t>& getPositionToAlignIds(const pos_t pos) const
    {
        return _positionToAlignIds[pos % MaxBufferSize];
    }
};
