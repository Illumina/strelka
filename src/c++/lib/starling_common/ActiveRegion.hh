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

#include "ActiveRegionReadBuffer.hh"
#include "CandidateSnvBuffer.hh"
#include "IndelBuffer.hh"
#include "alignment/GlobalAligner.hh"
#include "blt_util/align_path.hh"
#include "blt_util/blt_types.hh"
#include "options/IterativeAssemblerOptions.hh"
#include "starling_common/starling_types.hh"

#include <map>
#include <set>
#include <string>

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

    // Minimum supporting read count required to consider a haplotype for confirmation
    static const unsigned MinHaplotypeCount = 3u;

    /// Creates an active region object
    /// \param posRange position range of the active region
    /// \param ref reference
    /// \param maxIndelSize max indel size
    /// \param sampleIndex sample index
    /// \param aligner aligner for aligning haplotypes to the reference
    /// \param readBuffer read buffer
    /// \param indelBuffer indel buffer
    /// \param candidateSnvBuffer candidate SNV buffer
    /// \return active region object
    ActiveRegion(const pos_range& posRange,
                 const reference_contig_segment& ref,
                 const unsigned maxIndelSize,
                 const unsigned sampleIndex,
                 const GlobalAligner<int>& aligner,
                 const ActiveRegionReadBuffer& readBuffer,
                 IndelBuffer& indelBuffer,
                 CandidateSnvBuffer& candidateSnvBuffer):
        _posRange(posRange), _ref(ref), _maxIndelSize(maxIndelSize), _sampleIndex(sampleIndex),
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
    const pos_range _posRange;
    const reference_contig_segment& _ref;
    const unsigned _maxIndelSize;
    const unsigned _sampleIndex;
    const GlobalAligner<int> _aligner;

    const ActiveRegionReadBuffer& _readBuffer;
    IndelBuffer& _indelBuffer;
    CandidateSnvBuffer& _candidateSnvBuffer;

    std::set<align_id_t> _alignIdSoftClipped;

    /// Select the top haplotypes and convert these into primitive alleles
    ///
    /// \param[in] totalNumHaplotypingReads Total number of reads eligible for the haplotype generation process
    void
    processSelectedHaplotypes(
        HaplotypeToAlignIdSet& haplotypeToAlignIdSet,
        const unsigned totalNumHaplotypingReads);

    /// Create haplotypes using counting and process variants
    /// \return true if haplotype generation succeeds, false otherwise
    bool processHaplotypesWithCounting();

    /// Create haplotypes using assembly and process variants
    /// \return true if haplotype generation succeeds, false otherwise
    bool processHaplotypesWithAssembly();

    /// Do not use haplotyping to determine indel candidacy and MMDF relax positions
    void doNotUseHaplotyping();

    /// Convert the haplotype into primitive alleles and update _indelBuffer and _candidateSnvBuffer
    ///
    /// \param[in] totalNumHaplotypingReads Total number of reads eligible for the haplotype generation process, such
    ///               that alignIdList.size()/totalNumHaplotypingReads gives a meaning read count support ratio
    void convertToPrimitiveAlleles(
        const std::string& haploptypeSeq,
        const std::vector<align_id_t>& alignIdList,
        const unsigned totalNumHaplotypingReads,
        const uint8_t haplotypeId);
};
