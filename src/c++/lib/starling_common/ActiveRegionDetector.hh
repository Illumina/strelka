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

#include "ActiveRegionProcessor.hh"
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
/// Active regions are synchronized and shared by all samples. 
///
/// A major supporting component of the detector is the ActiveRegionReadBuffer, which does not store reads directly
/// but tracks variant and soft-clipping events per-position, and the associated read ids supporting each event.
///

class SampleActiveRegionDetector;

typedef known_pos_range2 ActiveRegion;

class ActiveRegionDetector
{
public:
    // alignment scores, same as bwa default values
    static const int ScoreMatch = 1;
    static const int ScoreMismatch = -4;
    static const int ScoreOpen = -5;
    static const int ScoreExtend = -1;
    static const int ScoreOffEdge = -100;

    /// Coordinates active region creation in all samples
    /// \param ref reference
    /// \param indelBuffer indel buffer
    /// \param candidateSnvBuffer candidate SNV buffer
    /// \param maxIndelSize maximum indel size
    /// \param sampleCount total sample count
    ActiveRegionDetector(
            const reference_contig_segment& ref,
            IndelBuffer& indelBuffer,
            CandidateSnvBuffer& candidateSnvBuffer,
            const unsigned maxIndelSize,
            const unsigned sampleCount);

    /// Gets an active region read buffer for the specified sample
    /// \param sampleIndex sample index
    ActiveRegionReadBuffer& getReadBuffer(const unsigned sampleIndex);

    /// Gets an active region ID for the specified position
    /// \param pos reference position
    /// \return active region ID or -1 if pos is not in an active region
    ActiveRegionId getActiveRegionId(const pos_t pos) const;

    /// Update the active region end position. Creates an active region if needed.
    /// \param pos reference position
    void updateEndPosition(const pos_t pos);

    /// Clear the detector. Create an active region if there exists an unclosed one.
    void clear();

    /// Clear the active region read buffer at pos
    void clearReadBuffer(const pos_t pos);

    /// Clear the position to active region map
    void clearPosToActiveRegionIdMapUpToPos(const pos_t pos);
private:
    const reference_contig_segment& _ref;
    const unsigned _sampleCount;
    std::vector<std::unique_ptr<SampleActiveRegionDetector>> _sampleActiveRegionDetector;
    IndelBuffer& _indelBuffer;
    CandidateSnvBuffer& _candidateSnvBuffer;
    const unsigned _maxIndelSize;
    // aligner to be used in active regions
    GlobalAligner<int> _aligner;

    ActiveRegion _synchronizedActiveRegion;
    RangeMap<pos_t, ActiveRegionId> _posToActiveRegionIdMap;

    SampleActiveRegionDetector& getSampleActiveRegionDetector(unsigned sampleIndex);
    void setPosToActiveRegionIdMap(ActiveRegion activeRegionRange);
    void updateActiveRegionRange(std::unique_ptr<ActiveRegion> activeRegion);
    void processExistingActiveRegion(const pos_t pos);
    void closeActiveRegion();
};

/// This class is responsible for detection of active regions for the specified sample
/// All the methods are accessed by ActiveRegionDetector
class SampleActiveRegionDetector
{
    friend class ActiveRegionDetector;

private:
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
        clearCoordinates();
    }

    /// Update the active region end position. Creates an active region if needed.
    /// \param pos reference position
    /// \return true if a new active region is created
    std::unique_ptr<ActiveRegion> updateEndPosition(const pos_t pos);

    /// Clear active region detector
    /// It may create an active region if there's unclosed one.
    /// \return true if a new active region is created
    std::unique_ptr<ActiveRegion> clear();

    /// Gets the read buffer
    /// \return read buffer
    ActiveRegionReadBuffer& getReadBuffer()
    {
        return _readBuffer;
    }

    void clearCoordinates();

    void clearReadBuffer(const pos_t pos);


    pos_t getActiveRegionStartPos() const
    {
        return _activeRegionStartPos;
    }

    std::unique_ptr<ActiveRegion> createActiveRegion();

    const reference_contig_segment& _ref;
    ActiveRegionReadBuffer _readBuffer;
    const unsigned _sampleIndex;

    /// true if updateEndPosition(pos) has not been called.
    bool _isBeginning;

    // An active region is created between [_activeRegionStartPos, _anchorPosFollowingPrevVariant
    pos_t _activeRegionStartPos;
    pos_t _anchorPosFollowingPrevVariant;

    /// Previous position classified as an anchor (anchors are non-STR loci that can be used to 'anchor' the edges
    /// of an active region)
    pos_t _prevAnchorPos;

    /// Previous position classified as a variant
    pos_t _prevVariantPos;

    /// The number of variants identified so far in the current candidate active region
    unsigned _numVariants;
};