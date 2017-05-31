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

#pragma once

#include "ActiveRegion.hh"
#include "blt_util/blt_types.hh"
#include "starling_read_segment.hh"
#include "indel.hh"
#include "IndelBuffer.hh"
#include "ActiveRegionReadBuffer.hh"
#include "CandidateSnvBuffer.hh"

#include <vector>
#include <list>
#include <set>


/// detects active regions
///
/// active regions are short genome segments where variation is sufficiently dense to trigger special haplotype handling methods
class ActiveRegionDetector
{
public:

    /// maximum buffer size in bases (must be larger than the maximum read size + max indel size
    static const unsigned MaxBufferSize = 1000u;

    /// max distance between two variants to be placed in the same active region
    static const unsigned MaxDistanceBetweenTwoVariants = 13u;

    /// min number of variants to form an active region
    static const unsigned MinNumVariantsPerRegion = 2u;

    // alignment scores, same as bwa default values
    static const int ScoreMatch = 1;
    static const int ScoreMismatch = -4;
    static const int ScoreOpen = -5;
    static const int ScoreExtend = -1;
    static const int ScoreOffEdge = -100;

    /// Creates an object that reads variant information and creates active regions
    /// \param ref reference segment
    /// \param indelBuffer indel buffer
    /// \param maxIndelSize maximum indel size
    /// \param sampleId sample Id
    ActiveRegionDetector(
        const reference_contig_segment& ref,
        IndelBuffer& indelBuffer,
        unsigned maxIndelSize,
        unsigned sampleId) :
        _ref(ref),
        _readBuffer(ref, indelBuffer),
        _indelBuffer(indelBuffer),
        _maxIndelSize(maxIndelSize),
        _sampleId(sampleId),
        _aligner(AlignmentScores<int>(ScoreMatch, ScoreMismatch, ScoreOpen, ScoreExtend, ScoreOffEdge, ScoreOpen, true, true))
    {
        _isBeginning = true;
        _activeRegionStartPos = -1;
        _anchorPosFollowingPrevVariant = 1;
        _prevAnchorPos = -1;
        _prevVariantPos = -1;
        _numVariants = 0;
    }

    /// Gets the read buffer
    /// \return read buffer
    ActiveRegionReadBuffer& getReadBuffer()
    {
        return _readBuffer;
    }

    const CandidateSnvBuffer& getCandidateSnvBuffer() const
    {
        return _candidateSnvBuffer;
    }

    ActiveRegionId getActiveRegionId(const pos_t pos) const
    {
        static const ActiveRegionId defaultActiveRegionId(-1);
        return _posToActiveRegionIdMap.getConstRefDefault(pos, defaultActiveRegionId);
    }

    void clearReadBuffer(const pos_t pos);

    /// update the active region end position. Creates an active region if needed.
    /// \param pos reference position
    void updateEndPosition(const pos_t pos);

    /// clear active region detector
    void clear();

    void clearUpToPos(const pos_t pos);

private:
    const reference_contig_segment& _ref;
    ActiveRegionReadBuffer _readBuffer;

    IndelBuffer& _indelBuffer;
    CandidateSnvBuffer _candidateSnvBuffer;

    const unsigned _maxIndelSize;
    const unsigned _sampleId;

    bool _isBeginning;
    pos_t _activeRegionStartPos;
    pos_t _anchorPosFollowingPrevVariant;
    pos_t _prevAnchorPos;
    pos_t _prevVariantPos;
    unsigned _numVariants;

    std::list<ActiveRegion> _activeRegions;

    // aligner to be used in active regions
    GlobalAligner<int> _aligner;

    RangeMap<pos_t, ActiveRegionId> _posToActiveRegionIdMap;

    void setPosToActiveRegionIdMap(pos_range activeRegionRange);
    void processActiveRegion();
};



