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

class SampleActiveRegionDetector
{
public:
    /// max distance between two variants to be placed in the same active region
    static const unsigned MaxDistanceBetweenTwoVariants = 13u;

    /// min number of variants to form an active region
    static const unsigned MinNumVariantsPerRegion = 2u;

    /// Creates an object that reads variant information and creates active regions
    /// \param ref reference segment
    /// \param indelBuffer indel buffer
    /// \param maxIndelSize maximum indel size
    /// \param sampleIndex sample Id
    SampleActiveRegionDetector(
        const reference_contig_segment& ref,
        IndelBuffer& indelBuffer,
        const unsigned sampleIndex) :
        _ref(ref),
        _readBuffer(ref, indelBuffer),
        _sampleIndex(sampleIndex)
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

    void clearReadBuffer(const pos_t pos);

    /// update the active region end position. Creates an active region if needed.
    /// \param pos reference position
    /// \return true if a new active region is created
    bool updateEndPosition(const pos_t pos);

    /// clear active region detector
    bool clear();

    known_pos_range2 getActiveRegionRange() const
    {
        return _activeRegionRange;
    }

    pos_t getActiveRegionStartPos() const
    {
        return _activeRegionStartPos;
    }

private:
    const reference_contig_segment& _ref;
    ActiveRegionReadBuffer _readBuffer;

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

    known_pos_range2 _activeRegionRange;

    void closeExistingActiveRegion();
};

class ActiveRegionDetector
{
public:
    // alignment scores, same as bwa default values
    static const int ScoreMatch = 1;
    static const int ScoreMismatch = -4;
    static const int ScoreOpen = -5;
    static const int ScoreExtend = -1;
    static const int ScoreOffEdge = -100;

    ActiveRegionDetector(
            const reference_contig_segment& ref,
            IndelBuffer& indelBuffer,
            CandidateSnvBuffer& candidateSnvBuffer,
            const unsigned maxIndelSize,
            const unsigned sampleCount) :
            _ref(ref),
            _sampleCount(sampleCount),
            _sampleActiveRegionDetector(sampleCount),
            _indelBuffer(indelBuffer),
            _candidateSnvBuffer(candidateSnvBuffer),
            _maxIndelSize(maxIndelSize),
            _aligner(AlignmentScores<int>(ScoreMatch, ScoreMismatch, ScoreOpen, ScoreExtend, ScoreOffEdge, ScoreOpen, true, true))

    {
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            _sampleActiveRegionDetector[sampleIndex].reset(
                    new SampleActiveRegionDetector(ref, indelBuffer, sampleIndex));
        }
    }

    ActiveRegionReadBuffer& getReadBuffer(const unsigned sampleIndex)
    {
        assert (sampleIndex < _sampleCount);
        return _sampleActiveRegionDetector[sampleIndex]->getReadBuffer();
    }

    ActiveRegionId getActiveRegionId(const pos_t pos) const
    {
        static const ActiveRegionId defaultActiveRegionId(-1);
        return _posToActiveRegionIdMap.getConstRefDefault(pos, defaultActiveRegionId);
    }

    void clearReadBuffer(const pos_t pos)
    {
        for (unsigned sampleIndex(0); sampleIndex<_sampleCount; ++sampleIndex)
        {
            _sampleActiveRegionDetector[sampleIndex]->clearReadBuffer(pos);
        }
    }

    /// update the active region end position. Creates an active region if needed.
    /// \param pos reference position
    void updateEndPosition(const pos_t pos)
    {
        processExistingActiveRegion(pos);
        for (unsigned sampleIndex(0); sampleIndex<_sampleCount; ++sampleIndex)
        {
            if (_sampleActiveRegionDetector[sampleIndex]->updateEndPosition(pos))
            {
                updateActiveRegionRange(sampleIndex);
            }
        }
    }

    /// clear active region detector
    void clear()
    {
        for (unsigned sampleIndex(0); sampleIndex<_sampleCount; ++sampleIndex)
        {
            if (_sampleActiveRegionDetector[sampleIndex]->clear())
            {
                updateActiveRegionRange(sampleIndex);
            }
        }
        if (_activeRegionRange.end_pos())
            closeActiveRegion();
    }

    void clearUpToPos(const pos_t pos);

    void updateActiveRegionRange(const unsigned sampleIndex)
    {
        const auto& sampleActiveRegionRange(getSampleActiveRegionDetector(sampleIndex).getActiveRegionRange());

        if (_activeRegionRange.end_pos())
        {
            // _activeRegionRange is valid
            _activeRegionRange.merge_range(sampleActiveRegionRange);
        }
        else
        {
            _activeRegionRange = sampleActiveRegionRange;
        }
    }

    void processExistingActiveRegion(const pos_t pos)
    {
        auto activeRegionEndPos(_activeRegionRange.end_pos());
        if (activeRegionEndPos &&
                (activeRegionEndPos <= (pos_t)(pos-ActiveRegionReadBuffer::MaxAssemblyPadding)))
        {
            bool canCloseActiveRegion(true);
            for (unsigned sampleIndex(0); sampleIndex<_sampleCount; ++sampleIndex)
            {
                auto unclosedActiveRegionStartPos(getSampleActiveRegionDetector(sampleIndex).getActiveRegionStartPos());
                if (unclosedActiveRegionStartPos >= 0 && unclosedActiveRegionStartPos < activeRegionEndPos-1)
                {
                    canCloseActiveRegion = false;
                    break;
                }
            }

            if (canCloseActiveRegion)
                closeActiveRegion();
        }
    }

    void closeActiveRegion()
    {
        assert (_activeRegionRange.size() > 0);
        for (unsigned sampleIndex(0); sampleIndex<_sampleCount; ++sampleIndex)
        {
            ActiveRegion activeRegion(_activeRegionRange,
                                      _ref, _maxIndelSize, sampleIndex,
                                      _aligner, _sampleActiveRegionDetector[sampleIndex]->getReadBuffer(),
                                      _indelBuffer, _candidateSnvBuffer);
            activeRegion.processHaplotypes();
        }
        setPosToActiveRegionIdMap(_activeRegionRange);
        _activeRegionRange.clear();
    }

private:
    const reference_contig_segment& _ref;
    const unsigned _sampleCount;
    std::vector<std::unique_ptr<SampleActiveRegionDetector>> _sampleActiveRegionDetector;
    IndelBuffer& _indelBuffer;
    CandidateSnvBuffer& _candidateSnvBuffer;
    const unsigned _maxIndelSize;
    // aligner to be used in active regions
    GlobalAligner<int> _aligner;

    known_pos_range2 _activeRegionRange;
    RangeMap<pos_t, ActiveRegionId> _posToActiveRegionIdMap;

    SampleActiveRegionDetector& getSampleActiveRegionDetector(unsigned sampleIndex)
    {
        assert(sampleIndex < _sampleCount);
        return *(_sampleActiveRegionDetector[sampleIndex]);
    }

    void setPosToActiveRegionIdMap(known_pos_range2 activeRegionRange);
};