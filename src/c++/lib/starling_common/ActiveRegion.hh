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

#include "blt_util/blt_types.hh"
#include "starling_common/starling_types.hh"
#include "alignment/GlobalAligner.hh"
#include "blt_util/align_path.hh"
#include "IndelBuffer.hh"
#include "ActiveRegionReadBuffer.hh"
#include "CandidateSnvBuffer.hh"

#include <string>
#include <map>
#include <set>
#include <options/SmallAssemblerOptions.hh>
#include <options/IterativeAssemblerOptions.hh>

typedef std::vector<RangeMap<pos_t,uint8_t>> RangeSet;
typedef std::map<std::string, std::vector<align_id_t>> HaplotypeToAlignIdSet;
typedef uint8_t HaplotypeId;

/// Represent all haplotypes found in the current active region
class ActiveRegion
{
public:
    // if the number of reads is larger than MinNumReadsToBypassAssembly, assembly is not conducted
    static const unsigned MinNumReadsToBypassAssembly = 1000u;

    // minimum fraction of reads covering the entire region to perform counting
    float MinFracReadsCoveringRegion = 0.65f;

    // if the region length is larger than MaxRefSpanToPerformAssembly, do not perform assembly
    static const unsigned MaxRefSpanToBypassAssembly = 250u;

    // assembly parameters
    const unsigned MaxAssemblyWordSize = 76u;
    const unsigned MinAssemblyCoverage = 3u;

    // minimum haplotype count to consider
    static const unsigned MinHaplotypeCount = 3u;

    /// Creates an active region object
    /// \param posRange position range of the active region
    /// \param ref reference
    /// \param maxIndelSize max indel size
    /// \param sampleCount sample count
    /// \param aligner aligner for aligning haplotypes to the reference
    /// \param readBuffer read buffer
    /// \param indelBuffer indel buffer
    /// \param candidateSnvBuffer candidate SNV buffer
    /// \return active region object
    ActiveRegion(const pos_range& posRange,
                 const reference_contig_segment& ref,
                 const unsigned maxIndelSize,
                 const unsigned sampleCount,
                 const GlobalAligner<int>& aligner,
                 const ActiveRegionReadBuffer& readBuffer,
                 IndelBuffer& indelBuffer,
                 CandidateSnvBuffer& candidateSnvBuffer):
        _posRange(posRange), _ref(ref), _maxIndelSize(maxIndelSize), _sampleCount(sampleCount),
        _aligner(aligner), _readBuffer(readBuffer), _indelBuffer(indelBuffer), _candidateSnvBuffer(candidateSnvBuffer)
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

    /// \param pos reference position
    /// \return true if pos belongs to this active region; false otherwise
    bool contains(const pos_t pos) const
    {
        return _posRange.is_pos_intersect(pos);
    }

    /// Decompose haplotypes into primitive alleles.
    /// Determine indel candidacy and register polymorphic sites to relax MMDF.
    void processHaplotypes();

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

    const ActiveRegionReadBuffer& _readBuffer;
    IndelBuffer& _indelBuffer;
    CandidateSnvBuffer& _candidateSnvBuffer;

    std::set<align_id_t> _alignIdSoftClipped;

    bool processSelectedHaplotypes(const unsigned sampleId, HaplotypeToAlignIdSet& haplotypeToAlignIdSet);

    /// Create haplotypes using counting and process variants
    /// \param sampleId sample id
    /// \return true if haplotype generation succeeds, false otherwise
    bool processHaplotypesWithCounting(const unsigned sampleId);

    /// Create haplotypes using assembly and process variants
    /// \param sampleId sample id
    /// \return true if haplotype generation succeeds, false otherwise
    bool processHaplotypesWithAssembly(const unsigned sampleId);

    /// Do not use haplotyping to determine indel candidacy and MMDF relax positions
    void doNotUseHaplotyping();

    /// convert the haplotype into primitive alleles and update _indelBuffer and _polySites
    void convertToPrimitiveAlleles(
        const unsigned sampleId,
        const std::string& haploptypeSeq,
        const std::vector<align_id_t>& alignIdList,
        const uint8_t haplotypeId);
};
