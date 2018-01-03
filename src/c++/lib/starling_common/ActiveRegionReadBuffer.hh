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

#pragma once

#include <vector>
#include "blt_util/reference_contig_segment.hh"
#include "blt_util/known_pos_range2.hh"
#include "starling_types.hh"
#include "indel.hh"
#include "IndelBuffer.hh"
#include "ReferenceRepeatFinder.hh"

/// AlignInfo object to store sample id and indel align type
struct AlignInfo
{
    unsigned sampleIndex;
    INDEL_ALIGN_TYPE::index_t indelAlignType;
    bool isForwardStrand;
};

struct ActiveRegionReadInfo
{
    /// Count of reads aligning to the active region, including reads which may not qualify for haplotype
    /// generation for various (configurable) reasons such as soft-clipping, partial coverage, etc...
    unsigned numReadsAlignedToActiveRegion;

    /// List of read segments which are eligible as haplotype generation input
    ///
    /// Invariant for this structure is (readSegmentsForHaplotypeGeneration.size() <= numReadsAlignedToActiveRegion)
    std::vector<std::pair<align_id_t, std::string>> readSegmentsForHaplotypeGeneration;
};

/// Helper object for ActiveRegionDetector: Tracks variant and anchor evidence per position, together with the read
/// ids supporting these
class ActiveRegionReadBuffer
{
public:

    /// Maximum buffer size in bases (must be larger than the maximum read size + max indel size
    static const unsigned MaxBufferSize = 1000u;

    /// Maximum read depth
    /// TODO: dynamically calculate maximum depth
    static const unsigned MaxDepth = 1000u;

    static const unsigned MinNumVariantsPerPosition = 9u;

    static const pos_t MaxAssemblyPadding = (pos_t)9u;

    /// Maximum repeat unit to consider
    static const unsigned MaxRepeatUnitLength = 50u;

    static const unsigned MinRepeatSpan = (pos_t)3u;

    /// Variant evidence weight for a mismatched basecall
    static const int MismatchWeight = 1;

    /// Variant evidence weight for an indel
    static const int IndelWeight = 4;

    // if the fraction is larger than MinAlternativeAlleleFractionLowDepth
    // the position becomes candidate even if the number is lower than MinNumVariantsPerPosition
    const float MinAlternativeAlleleFractionLowDepth = 0.35;

    /// Read buffer to be used in active regions
    /// \param ref reference
    /// \param indelBuffer indel buffer
    ActiveRegionReadBuffer(
        const reference_contig_segment& ref,
        const float minAlternativeAlleleFraction,
        IndelBuffer& indelBuffer)
        :
        _ref(ref),
        _minAlternativeAlleleFraction(minAlternativeAlleleFraction),
        _refRepeatFinder(ref, MaxRepeatUnitLength, MaxBufferSize, MinRepeatSpan),
        _indelBuffer(indelBuffer),
        _variantCounter(MaxBufferSize),
        _depth(MaxBufferSize),
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
    /// \param segmentSeq the soft-clip segment sequence
    /// \param isBeginEdge true if the soft-clip is at the beginning
    void insertSoftClipSegment(const align_id_t alignId, const pos_t pos, const std::string& segmentSeq, const bool isBeginEdge);

    /// insert indel
    /// \param indelObservation indel observation object
    void insertIndel(const IndelObservation& indelObservation);

    /// checks if pos is an anchor position
    /// \param pos reference position
    /// \return true if pos is an anchor position, false otherwise
    bool isAnchor(const pos_t pos) const
    {
        return _refRepeatFinder.isAnchor(pos);
    }

    /// Set end position of the buffer
    /// \param endPos end position (exclusive)
    void setEndPos(const pos_t endPos);

    /// Gets the beginning position
    /// \return begin position
    pos_t getBeginPos() const
    {
        return _readBufferRange.begin_pos();
    }

    /// Gets the end position
    /// \return end position
    pos_t getEndPos() const
    {
        return _readBufferRange.end_pos();
    }

    /// Gets sample id and indel align type
    /// \param alignId align id
    /// \return sample id and indel align type
    const AlignInfo& getAlignInfo(const align_id_t alignId) const
    {
        return _alignIdToAlignInfo[alignId % MaxDepth];
    }

    /// Gets read segments for the range
    /// \param posRange position range
    /// \param readInfo read info object to store read segments
    /// \param includePartialReads if true, only reads fully covering the region will be retrieved
    void getReadSegments(
        const known_pos_range2& posRange,
        ActiveRegionReadInfo& readInfo,
        const bool includePartialReads,
        const unsigned minReadSegmentLength = 1u) const;

    /// cache sampleIndex and indelAlignType corresponding to alignId
    /// \param alignId align id
    /// \param sampleIndex sample index
    /// \param indelAlignType indel align type
    void setAlignInfo(const align_id_t alignId, const unsigned sampleIndex, const INDEL_ALIGN_TYPE::index_t indelAlignType, const bool isForwardStrand)
    {
        AlignInfo& alignInfo = _alignIdToAlignInfo[alignId % MaxDepth];
        alignInfo.sampleIndex = sampleIndex;
        alignInfo.indelAlignType = indelAlignType;
        alignInfo.isForwardStrand = isForwardStrand;
    }

    /// Gets sample id
    /// \param alignId align id
    /// \return sample id
    unsigned getSampleIndex(const align_id_t alignId) const
    {
        return _alignIdToAlignInfo[alignId % MaxDepth].sampleIndex;
    }

    /// Clear buffer
    /// \param pos position
    void clearPos(const pos_t pos)
    {
        _positionToAlignIds[pos % MaxBufferSize].clear();
        resetCounter(pos);
        _readBufferRange.set_begin_pos(pos+1);
    }

    /// Checks if pos is candidate variant position
    /// \param pos position
    /// \return true if pos is candidate variant, false otherwise
    bool isCandidateVariant(const pos_t pos) const;

    /// Checks if the depth is zero at pos
    bool isDepthZero(const pos_t pos) const
    {
        return (getDepth(pos) == 0u);
    }

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
    const float _minAlternativeAlleleFraction;

    ReferenceRepeatFinder _refRepeatFinder;

    IndelBuffer& _indelBuffer;

    known_pos_range2 _readBufferRange;

    /// Stores the number of variant observations (including soft-clip) per position
    std::vector<unsigned> _variantCounter;

    /// Store the number of non-clipped observations per position
    std::vector<unsigned> _depth;

    // for haplotypes
    std::vector<std::vector<align_id_t>> _positionToAlignIds;

    // to store align information
    std::vector<AlignInfo> _alignIdToAlignInfo;

    std::vector<std::vector<VariantType>> _variantInfo;
    std::vector<std::vector<std::string>> _insertSeqBuffer;
    char _snvBuffer[MaxDepth][MaxBufferSize];

    void setMatch(const align_id_t id, const pos_t pos);
    void setMismatch(const align_id_t id, const pos_t pos, const char baseChar);
    void setDelete(const align_id_t id, const pos_t pos);
    void setInsert(const align_id_t id, const pos_t pos, const std::string& insertSeq);
    void setSoftClipSegment(const align_id_t id, const pos_t pos, const std::string& segmentSeq);
    bool getHaplotypeBase(const align_id_t id, const pos_t pos, std::string& base) const;

    void resetCounter(const pos_t pos)
    {
        int index = pos % MaxBufferSize;

        _variantCounter[index] = 0;
        _depth[index] = 0;
    }

    void addVariantCount(const pos_t pos, const unsigned count)
    {
        int index = pos % MaxBufferSize;
        _variantCounter[index] += count;
        ++_depth[index];
    }

    void addSoftClipCount(const pos_t pos, const unsigned count)
    {
        int index = pos % MaxBufferSize;
        _variantCounter[index] += count;
    }

    unsigned getVariantCount(const pos_t pos) const
    {
        return _variantCounter[pos % MaxBufferSize];
    }

    void addAlignIdToPos(const align_id_t alignId, const pos_t pos)
    {
        int index = pos % MaxBufferSize;
        if (_positionToAlignIds[index].empty() || _positionToAlignIds[index].back() != alignId)
            _positionToAlignIds[index].push_back(alignId);
    }

    unsigned getDepth(const pos_t pos) const
    {
        int index = pos % MaxBufferSize;
        return _depth[index];
    }

    const std::vector<align_id_t>& getPositionToAlignIds(const pos_t pos) const
    {
        return _positionToAlignIds[pos % MaxBufferSize];
    }
};
