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
#include "ActiveRegionReadBuffer.hh"

#include <vector>
#include <list>
#include <set>

/// detects active regions
///
/// active regions are short genome segmetns where variation is sufficiently dense to trigger special haplotype handling methods
class ActiveRegionDetector
{
public:

    // maximum buffer size in bases (must be larger than the maximum read size + max indel size
    static const unsigned MaxBufferSize = 1000u;

    // minimum read depth
    static const unsigned MinDepth = ActiveRegion::MinHaplotypeCount;


    // maximum distance between two variants belonging to the same active region
    static const int MaxDistanceBetweenTwoVariants = 13u;

    // alignment scores, same as bwa default values
    static const int ScoreMatch = 1;
    static const int ScoreMismatch = -4;
    static const int ScoreOpen = -5;
    static const int ScoreExtend = -1;
    static const int ScoreOffEdge = -100;

    // for expansion of active regions
    const pos_t MaxRepeatSpan = 100u;

    // minimum alternative allele fraction to call a position as a candidate variant
    const float MinAlternativeAlleleFraction = 0.2;

    /// Creates an object that reads variant information and creates active regions
    /// \param ref reference segment
    /// \param indelBuffer indel buffer
    /// \param maxIndelSize maximum indel size
    /// \param sampleCount sample count
    /// \param maxDetectionWindowSize maximum active region size
    /// \param minNumVariantsPerPosition minimum number of variants per position to make the position as a candidate variant pos
    /// \param minNumVariantsPerRegion minimum number of variants per region to create an active region
    ActiveRegionDetector(
        const reference_contig_segment& ref,
        IndelBuffer& indelBuffer,
        unsigned maxIndelSize,
        unsigned sampleCount,
        unsigned maxDetectionWindowSize = 100,
        unsigned minNumVariantsPerPosition = 9,
        unsigned minNumVariantsPerRegion = 2) :
        _ref(ref),
        _readBuffer(ref, sampleCount, indelBuffer),
        _indelBuffer(indelBuffer),
        _maxIndelSize(maxIndelSize),
        _sampleCount(sampleCount),
        _maxDetectionWindowSize(maxDetectionWindowSize),
        _minNumVariantsPerPositionPerSample(minNumVariantsPerPosition),
        _minNumVariantsPerRegion(minNumVariantsPerRegion),
        _polySites(sampleCount),
        _aligner(AlignmentScores<int>(ScoreMatch, ScoreMismatch, ScoreOpen, ScoreExtend, ScoreOffEdge, ScoreOpen, true, true)),
        _alignerForAssembly(AlignmentScores<int>(ScoreMatch, ScoreMismatch, ScoreOpen, ScoreExtend, 0, ScoreOpen, true, true))
    {
        _numVariants = 0;
        _activeRegionPtr = nullptr;
        _activeRegionStartPos = 0;
        _prevVariantPos = 0;
        _lastActiveRegionEnd = 0;
    }

    ActiveRegionReadBuffer& getReadBuffer()
    {
        return _readBuffer;
    }

    /// update the active region buffer start position
    /// \param pos reference position
    void updateStartPosition(const pos_t pos);

    /// update the active region end position. Creates an active region if needed.
    /// \param pos reference position
    /// \param isLastPos
    void updateEndPosition(const pos_t pos);

    /// Checks if mismatches occur consistently at position pos
    /// \param sampleId sample id
    /// \param pos reference position
    /// \return true if pos is a polymorphic site; false otherwise.
    bool isPolymorphicSite(const unsigned sampleId, const pos_t pos) const;

    void clear();

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
    ActiveRegionReadBuffer _readBuffer;

    IndelBuffer& _indelBuffer;

    const unsigned _maxIndelSize;
    const unsigned _sampleCount;
    const unsigned _maxDetectionWindowSize;

    const unsigned _minNumVariantsPerPositionPerSample;
    const unsigned _minNumVariantsPerRegion;

    pos_t _activeRegionStartPos;
    pos_t _prevVariantPos;
    pos_t _lastActiveRegionEnd;
    unsigned _numVariants;

    std::unique_ptr<ActiveRegion> _activeRegionPtr;

    // record polymorphic sites
    RangeSet _polySites;

    // aligner to be used in active regions
    GlobalAligner<int> _aligner;
    GlobalAligner<int> _alignerForAssembly;

    // expand the active region to cover repeats
    void getExpandedRange(const pos_range& origActiveRegion, pos_range& newActiveRegion);

    void processActiveRegion();

    bool isCandidateVariant(const pos_t pos) const;
    bool isInvariant(const pos_t pos) const;
};



