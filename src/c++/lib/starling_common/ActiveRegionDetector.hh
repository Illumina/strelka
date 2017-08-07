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

#pragma once

#include "ActiveRegion.hh"
#include "blt_util/blt_types.hh"
#include "starling_read_segment.hh"
#include "IndelBuffer.hh"
#include "ActiveRegionReadBuffer.hh"
#include "CandidateSnvBuffer.hh"

#include <list>


/// \brief Detects active regions
///
/// Active regions are short genome segments where variation is sufficiently dense to trigger special haplotype
/// detection and handling methods
///
/// Active regions are defined so as to cluster a group of nearby variant candidates so that they can be handled by
/// haplotype analysis. Given a particular cluster of variant the active region is extended so that it starts and
/// ends at an 'anchor' point, a non-STR region appropriate for defining the endpoint of an active region so as to
/// avoid common artifacts associated with partial representation of an STR.
///
/// A major supporting component of the detector is the ActiveRegionReadBuffer, which does not store reads directly
/// but tracks variant and soft-clipping events per-position, and the associated read ids supporting each event.
///
class ActiveRegionDetector
{
public:
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
    /// \param sampleIndex sample Id
    ActiveRegionDetector(
        const reference_contig_segment& ref,
        IndelBuffer& indelBuffer,
        CandidateSnvBuffer& candidateSnvBuffer,
        unsigned maxIndelSize,
        unsigned sampleIndex) :
        _ref(ref),
        _readBuffer(ref, indelBuffer),
        _indelBuffer(indelBuffer),
        _candidateSnvBuffer(candidateSnvBuffer),
        _maxIndelSize(maxIndelSize),
        _sampleIndex(sampleIndex),
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
    CandidateSnvBuffer& _candidateSnvBuffer;

    const unsigned _maxIndelSize;
    const unsigned _sampleIndex;

    bool _isBeginning;
    pos_t _activeRegionStartPos;
    pos_t _anchorPosFollowingPrevVariant;

    /// Previous position classified as an anchor (anchors are non-STR loci that can be used to 'anchor' the edges
    /// of an active region)
    pos_t _prevAnchorPos;

    /// Previous position classified as a variant
    pos_t _prevVariantPos;

    /// The number of variants identified so far in the current candidate active region
    unsigned _numVariants;

    /// \TODO Why does the object support multiple active regions when processActiveRegions will only process one. Why is this a list? (STREL-655)
    std::list<ActiveRegion> _activeRegions;

    // aligner to be used in active regions
    GlobalAligner<int> _aligner;

    RangeMap<pos_t, ActiveRegionId> _posToActiveRegionIdMap;

    void setPosToActiveRegionIdMap(pos_range activeRegionRange);
    void processActiveRegion();

    void closeExistingActiveRegion();
};

