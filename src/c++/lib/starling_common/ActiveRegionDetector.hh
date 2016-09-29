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

#include "ActiveRegion.hh"
#include "blt_util/blt_types.hh"
#include "starling_read_segment.hh"
#include "indel.hh"
#include "IndelBuffer.hh"

#include <vector>
#include <list>
#include <set>

/// An agent that detects active regions
class ActiveRegionDetector
{
public:

    // maximum buffer size in bases (must be larger than the maximum read size + max indel size
    static const unsigned MaxBufferSize = 1000u;

    // maximum read depth
    static const unsigned MaxDepth = ActiveRegion::MaxDepth;

    // minimum read depth
    static const unsigned MinDepth = ActiveRegion::MinHaplotypeCount;

    // variant count to add for a single mismatch or indel
    static const int MismatchWeight = 1;
    static const int IndelWeight = 4;

    // maximum distance between two variants belonging to the same active region
    static const int MaxDistanceBetweenTwoVariants = 14u;

    // alignment scores, same as bwa default values
    static const int ScoreMatch = 1;
    static const int ScoreMismatch = -4;
    static const int ScoreOpen = -5;
    static const int ScoreExtend = -1;
    static const int ScoreOffEdge = -100;

    // for expansion of active regions
    const unsigned MaxRepeatUnitLength = 3u;
    const pos_t MaxRepeatSpan = 20u;

    // minimum alternative allele fraction to call a position as a candidate variant
    const float MinAlternativeAlleleFraction = 0.2;

    /// Creates an object that reads variant information and creates active regions
    /// \param ref reference segment
    /// \param indelBuffer indel buffer
    /// \param maxDeletionSize maximum deletion size
    /// \param sampleCount sample count
    /// \param maxDetectionWindowSize maximum active region size
    /// \param minNumVariantsPerPosition minimum number of variants per position to make the position as a candidate variant pos
    /// \param minNumVariantsPerRegion minimum number of variants per region to create an active region
    /// \return active region detector object
    ActiveRegionDetector(
        const reference_contig_segment& ref,
        IndelBuffer& indelBuffer,
        unsigned maxDeletionSize,
        unsigned sampleCount,
        unsigned maxDetectionWindowSize = 30,
        unsigned minNumVariantsPerPosition = 9,
        unsigned minNumVariantsPerRegion = 2) :
        _ref(ref),
        _indelBuffer(indelBuffer),
        _maxDeletionSize(maxDeletionSize),
        _sampleCount(sampleCount),
        _maxDetectionWindowSize(maxDetectionWindowSize),
        _minNumVariantsPerPositionPerSample(minNumVariantsPerPosition),
        _minNumVariantsPerRegion(minNumVariantsPerRegion),
        _variantCounter(sampleCount, std::vector<unsigned>(MaxBufferSize)),
        _depth(sampleCount, std::vector<unsigned>(MaxBufferSize)),
        _positionToAlignIds(MaxBufferSize),
        _alignIdToAlignInfo(MaxDepth),
        _variantInfo(MaxDepth, std::vector<VariantType>(MaxBufferSize, VariantType())),
        _insertSeqBuffer(MaxDepth, std::vector<std::string>(MaxBufferSize, std::string())),
        _polySites(sampleCount),
        _aligner(AlignmentScores<int>(ScoreMatch, ScoreMismatch, ScoreOpen, ScoreExtend, ScoreOffEdge, ScoreOpen, true, true))
    {
        _numVariants = 0;
        _activeRegionStartPos = 0;
        _prevVariantPos = 0;
        _lastActiveRegionEnd = 0;
    }

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
    void insertSoftClipSegment(const align_id_t alignId, const pos_t pos, const std::string& segmentSeq);

    /// insert indel
    /// \param sampleId sample id
    /// \param indelObservation indel observation object
    void insertIndel(const unsigned sampleId, const IndelObservation& indelObservation);

    /// update the active region buffer start position
    /// \param pos reference position
    void updateStartPosition(const pos_t pos);

    /// update the active region end position. Creates an active region if needed.
    /// \param pos reference position
    /// \param isLastPos
    void updateEndPosition(const pos_t pos, const bool isLastPos);

    /// cache sampleId and indelAlignType corresponding to alignId
    /// \param alignId align id
    /// \param sampleId sample id
    /// \param indelAlignType indel align type
    inline void setAlignInfo(const align_id_t alignId, unsigned sampleId, INDEL_ALIGN_TYPE::index_t indelAlignType)
    {
        AlignInfo& alignInfo = _alignIdToAlignInfo[alignId % MaxDepth];
        alignInfo.sampleId = sampleId;
        alignInfo.indelAlignType = indelAlignType;
    }

    inline unsigned getSampleId(const align_id_t alignId)
    {
        return _alignIdToAlignInfo[alignId % MaxDepth].sampleId;
    }

    /// Checks if mismatches occur consistently at position pos
    /// \param sampleId sample id
    /// \param pos reference position
    /// \return true if pos is a polymorphic site; false otherwise.
    bool isPolymorphicSite(const unsigned sampleId, const pos_t pos) const;

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
    IndelBuffer& _indelBuffer;

    const unsigned _maxDeletionSize;
    const unsigned _sampleCount;
    const unsigned _maxDetectionWindowSize;

    const unsigned _minNumVariantsPerPositionPerSample;
    const unsigned _minNumVariantsPerRegion;

    pos_t _activeRegionStartPos;
    pos_t _prevVariantPos;
    pos_t _lastActiveRegionEnd;
    unsigned _numVariants;

    std::list<ActiveRegion> _activeRegions;
    std::vector<std::vector<unsigned>> _variantCounter;
    std::vector<std::vector<unsigned>> _depth;

    // for haplotypes
    std::vector<std::vector<align_id_t>> _positionToAlignIds;

    // to store align information
    std::vector<AlignInfo> _alignIdToAlignInfo;

    std::vector<std::vector<VariantType>> _variantInfo;
    std::vector<std::vector<std::string>> _insertSeqBuffer;
    char _snvBuffer[MaxDepth][MaxBufferSize];

    // record polymorphic sites
    RangeSet _polySites;

    // aligner to be used in active regions
    GlobalAligner<int> _aligner;

    // expand the active region to cover repeats
    void getExpandedRange(const pos_range& origActiveRegion, pos_range& newActiveRegion);

    bool isCandidateVariant(const pos_t pos) const;

    inline void resetCounter(const pos_t pos)
    {
        int index = pos % MaxBufferSize;

        for (unsigned sampleId=0; sampleId<_sampleCount; ++sampleId)
        {
            _variantCounter[sampleId][index] = 0;
            _depth[sampleId][index] = 0;
        }
    }

    inline void addVariantCount(const unsigned sampleId, const pos_t pos, unsigned count)
    {
        int index = pos % MaxBufferSize;
        _variantCounter[sampleId][index] += count;
        ++_depth[sampleId][index];
    }

    inline unsigned getVariantCount(const unsigned sampleId, const pos_t pos) const
    {
        return _variantCounter[sampleId][pos % MaxBufferSize];
    }

    inline void addAlignIdToPos(const align_id_t alignId, const pos_t pos)
    {
        int index = pos % MaxBufferSize;
        if (_positionToAlignIds[index].empty() || _positionToAlignIds[index].back() != alignId)
            _positionToAlignIds[index].push_back(alignId);
    }

    inline unsigned getDepth(const unsigned sampleId, const pos_t pos) const
    {
        int index = pos % MaxBufferSize;
        return _depth[sampleId][index];
    }

    inline const std::vector<align_id_t>& getPositionToAlignIds(const pos_t pos) const
    {
        return _positionToAlignIds[pos % MaxBufferSize];
    }

    void setMatch(const align_id_t id, const pos_t pos);
    void setMismatch(const align_id_t id, const pos_t pos, char baseChar);
    void setDelete(const align_id_t id, const pos_t pos);
    void setInsert(const align_id_t id, const pos_t pos, const std::string& insertSeq);
    void setSoftClipSegment(const align_id_t id, const pos_t pos, const std::string& segmentSeq);

    bool setHaplotypeBase(const align_id_t id, const pos_t pos, std::string& base) const;

    inline void clearPos(pos_t pos)
    {
        _positionToAlignIds[pos % MaxBufferSize].clear();
        resetCounter(pos);
    }
};



