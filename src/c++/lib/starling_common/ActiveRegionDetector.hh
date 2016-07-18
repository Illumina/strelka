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
    static const unsigned MaxBufferSize = 1000;
    static const unsigned MaxDepth = 1000;

    // alignment scores, same as bwa default values
    static const int ScoreMatch = 1;
    static const int ScoreMismatch = -4;
    static const int ScoreOpen = -5;
    static const int ScoreExtend = -1;
    static const int ScoreOffEdge = -100;

    ActiveRegionDetector(
        const reference_contig_segment& ref,
        IndelBuffer& indelBuffer,
        unsigned maxDeletionSize,
        unsigned maxDetectionWindowSize = 30,
        unsigned minNumMismatchesPerPosition = 9,
        unsigned minNumVariantsPerRegion = 2) :
        _ref(ref),
        _indelBuffer(indelBuffer),
        _maxDeletionSize(maxDeletionSize),
        _maxDetectionWindowSize(maxDetectionWindowSize),
        _minNumMismatchesPerPosition(minNumMismatchesPerPosition),
        _minNumVariantsPerRegion(minNumVariantsPerRegion),
        _variantCounter(MaxBufferSize),
        _positionToAlignIds(MaxBufferSize),
        _variantInfo(MaxDepth, std::vector<VariantType>(MaxBufferSize, VariantType())),
        _insertSeqBuffer(MaxDepth, std::vector<std::string>(MaxBufferSize, std::string())),
        _aligner(AlignmentScores<int>(ScoreMatch, ScoreMismatch, ScoreOpen, ScoreExtend, ScoreOffEdge, ScoreOpen, true, true))
    {
        _bufferStartPos = 0;

        _numVariants = 0;
        _activeRegionStartPos = 0;
        _prevVariantPos = 0;
    }

    void insertMatch(const align_id_t alignId, const pos_t pos);
    void insertMismatch(const align_id_t alignId, const pos_t pos, const char baseChar);
    void insertIndel(const IndelObservation& indelObservation);
    void updateStartPosition(const pos_t pos);
    void updateEndPosition(const pos_t pos);
    bool isEmpty() const
    {
        return _activeRegions.empty();
    }
    bool isPolymorphicSite(const pos_t pos) const;

    enum VariantType
    {
        MATCH,
        MISMATCH,
        DELETE,
        INSERT,
        MISMATCH_INSERT
    };

private:
    const reference_contig_segment& _ref;
    IndelBuffer& _indelBuffer;

    unsigned _maxDeletionSize;
    unsigned _maxDetectionWindowSize;
    unsigned _minNumMismatchesPerPosition;
    unsigned _minNumVariantsPerRegion;

    pos_t _bufferStartPos;

    pos_t _prevVariantPos;
    pos_t _activeRegionStartPos;
    unsigned _numVariants;

    std::list<ActiveRegion> _activeRegions;
    std::vector<unsigned> _variantCounter;

    // for haplotypes
    std::vector<std::vector<align_id_t>> _positionToAlignIds;
    std::vector<std::vector<VariantType>> _variantInfo;
    std::vector<std::vector<std::string>> _insertSeqBuffer;
    char _snvBuffer[MaxDepth][MaxBufferSize];

    // record polymorphic sites
    RangeSet _polySites;

    // aligner to be used in active regions
    GlobalAligner<int> _aligner;

    bool isCandidateVariant(const pos_t pos) const;

    inline void resetCounter(const pos_t pos)
    {
        _variantCounter[pos % MaxBufferSize] = 0;
    }

    inline void addCount(const pos_t pos, unsigned count = 1)
    {
        _variantCounter[pos % MaxBufferSize] += count;
    }

    inline unsigned getCount(const pos_t pos) const
    {
        return _variantCounter[pos % MaxBufferSize];
    }

    inline void addAlignIdToPos(const align_id_t alignId, const pos_t pos)
    {
        int index = pos % MaxBufferSize;
        if (_positionToAlignIds[index].empty() || _positionToAlignIds[index].back() != alignId)
            _positionToAlignIds[index].push_back(alignId);
    }

    inline const std::vector<align_id_t>& getPositionToAlignIds(const pos_t pos) const
    {
        return _positionToAlignIds[pos % MaxBufferSize];
    }

    void setMatch(const align_id_t id, const pos_t pos);
    void setMismatch(const align_id_t id, const pos_t pos, char baseChar);
    void setDelete(const align_id_t id, const pos_t pos);
    void setInsert(const align_id_t id, const pos_t pos, const std::string& insertSeq);
    void setHaplotypeBase(const align_id_t id, const pos_t pos, std::string& base) const;

    inline void clearPos(pos_t pos)
    {
        _positionToAlignIds[pos % MaxBufferSize].clear();
    }
};



