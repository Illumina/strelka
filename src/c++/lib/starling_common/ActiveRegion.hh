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

#include <string>
#include <map>
#include <set>

typedef RangeMap<pos_t,unsigned char> RangeSet;

/// AlignInfo object to store sample id and indel align type
struct AlignInfo
{
    AlignInfo() {}

    unsigned sampleId;
    INDEL_ALIGN_TYPE::index_t indelAlignType;
};

/// Represent all haplotypes found in the current active region
class ActiveRegion
{
public:
    // maximum read depth
    static const unsigned MaxDepth = 1000;

    // minimum haplotype count to consider
    static const unsigned MinHaplotypeCount = 3;

    // minimum haplotype frequency to relax MMDF
    const float HaplotypeFrequencyThreshold = 0.4;

    const char missingPrefix = '.';

    /// Creates an active region object
    /// \param start start position
    /// \param end end position
    /// \param refSeq reference sequence corresponding to this active region
    /// \param aligner aligner for aligning haplotypes to the reference
    /// \param alignIdToAlignInfo map from align id to (sampleId, indelAlignType)
    /// \return active region object
    ActiveRegion(pos_t start, pos_t end, const std::string& refSeq, const GlobalAligner<int>& aligner, const std::vector<AlignInfo>& alignIdToAlignInfo):
        _start(start), _end(end), _refSeq(refSeq),
        _aligner(aligner),
        _alignIdToAlignInfo(alignIdToAlignInfo),
        _alignIdToHaplotype(),
        _alignIdReachingEnd()
    {}

    /// \return start position of the active region
    pos_t getStart() const
    {
        return _start;
    }

    /// \param pos reference position
    /// \return true if pos belongs to this active region; false otherwise
    bool contains(pos_t pos) const
    {
        return pos >= _start && pos <= _end;
    }

    /// Insert haplotype base of the alignment into the position
    /// \param alignId align id
    /// \param pos reference position
    /// \param base base at position pos. base.length() is 1 for match/mismatch, >1 for insertion, and 0 for deletion.
    void insertHaplotypeBase(align_id_t alignId, pos_t pos, const std::string& base);

    /// Decompose haplotypes into primitive alleles.
    /// Determine indel candidacy and regiser polymorphic sites to relax MMDF.
    /// \param indelBuffer
    /// \param polySites
    void processHaplotypes(IndelBuffer& indelBuffer, RangeSet& polySites) const;

    void setSoftClipped(const align_id_t alignId)
    {
        _alignIdSoftClipped.insert(alignId);
    }

private:
    pos_t _start;
    pos_t _end;
    const std::string& _refSeq;
    const GlobalAligner<int> _aligner;
    const std::vector<AlignInfo>&  _alignIdToAlignInfo;

    std::map<align_id_t, std::string> _alignIdToHaplotype;
    std::set<align_id_t> _alignIdReachingEnd;
    std::set<align_id_t> _alignIdSoftClipped;

    void convertToPrimitiveAlleles(
        const std::string& haploptypeSeq,
        const std::vector<align_id_t>& alignIdList,
        const unsigned totalReadCount,
        IndelBuffer& indelBuffer,
        RangeSet& polySites) const;
};
