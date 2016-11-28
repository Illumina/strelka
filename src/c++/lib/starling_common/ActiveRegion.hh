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

#include "blt_util/blt_types.hh"
#include "starling_common/starling_types.hh"
#include "alignment/GlobalAligner.hh"
#include "blt_util/align_path.hh"
#include "IndelBuffer.hh"
#include "ActiveRegionReadBuffer.hh"

#include <string>
#include <map>
#include <set>
#include <options/SmallAssemblerOptions.hh>
#include <options/IterativeAssemblerOptions.hh>

typedef std::vector<RangeMap<pos_t,unsigned char>> RangeSet;

/// Represent all haplotypes found in the current active region
class ActiveRegion
{
public:
    // if haplotype depth is at least HighDepth, low depth filter is not applied
    static const unsigned HighDepth = 20u;

    // minimum haplotype count to consider
    static const unsigned MinHaplotypeCount = 3u;
    static const unsigned MaxSNVHpolSize = 4u;

    // minimum haplotype frequency to relax MMDF
    const float HaplotypeFrequencyThreshold = 0.4;
    const unsigned AssemblyFlankingSequenceLength = 10u;

    /// Creates an active region object
    /// \param posRange position range of the active region
    /// \param ref reference
    /// \param aligner aligner for aligning haplotypes to the reference
    /// \param readBuffer read buffer
    /// \return active region object
    ActiveRegion(pos_range& posRange,
                 const reference_contig_segment& ref,
                 const unsigned maxIndelSize,
                 const unsigned sampleCount,
                 const GlobalAligner<int>& aligner,
                 const GlobalAligner<int>& alignerForAssembly,
                 const ActiveRegionReadBuffer& readBuffer):
        _posRange(posRange), _ref(ref), _maxIndelSize(maxIndelSize), _sampleCount(sampleCount),
        _aligner(aligner), _alignerForAssembly(alignerForAssembly),
        _readBuffer(readBuffer)
    {
    }

    /// \return begin position of the active region
    pos_t getBeginPosition() const
    {
        return _posRange.begin_pos;
    }

    /// \return begin position of the active region
    pos_t getEndPosition() const
    {
        return _posRange.end_pos;
    }

    void extendEndPosition(const pos_t endPos)
    {
        _posRange.set_end_pos(endPos);
    }

    /// \param pos reference position
    /// \return true if pos belongs to this active region; false otherwise
    bool contains(pos_t pos) const
    {
        return _posRange.is_pos_intersect(pos);
    }

    /// Decompose haplotypes into primitive alleles.
    /// Determine indel candidacy and regiser polymorphic sites to relax MMDF.
    /// \param indelBuffer
    /// \param polySites
    void processHaplotypes(IndelBuffer &indelBuffer, RangeSet &polySites);

    /// Mark a read soft-clipped
    /// \param alignId align id
    void setSoftClipped(const align_id_t alignId)
    {
        _alignIdSoftClipped.insert(alignId);
    }

private:
    pos_range _posRange;
    const reference_contig_segment& _ref;
    const unsigned _maxIndelSize;
    const unsigned _sampleCount;
    const GlobalAligner<int> _aligner;
    const GlobalAligner<int> _alignerForAssembly;

    const ActiveRegionReadBuffer& _readBuffer;
    std::set<align_id_t> _alignIdReachingStart;
    std::set<align_id_t> _alignIdReachingEnd;
    std::set<align_id_t> _alignIdSoftClipped;

    bool processHaplotypesWithCounting(IndelBuffer& indelBuffer, RangeSet& polySites, unsigned sampleId) const;
    void processHaplotypesWithAssembly(IndelBuffer& indelBuffer, RangeSet& polySites, unsigned sampleId) const;

    void convertToPrimitiveAlleles(
        const unsigned sampleId,
        const std::string& haploptypeSeq,
        const std::vector<align_id_t>& alignIdList,
        const unsigned totalReadCount,
        const bool isTopTwo,
        const pos_range posRange,
        IndelBuffer& indelBuffer,
        RangeSet& polySites) const;
};
