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

/// A worker that
/// (1) creates haplotypes,
/// (2) align them to the reference to discover primitive alleles, and
/// (3) put indels/SNVs into IndelBuffer and CandidateSnvBuffer for downstream processing
class ActiveRegionProcessor
{
public:
    /// if the number of reads is larger than MinNumReadsToBypassAssembly, assembly is not conducted
    static const unsigned MinNumReadsToBypassAssembly = 1000u;

    /// minimum fraction of reads covering the entire region to perform counting
    float MinFracReadsCoveringRegion = 0.65f;

    /// if the region length is larger than MaxRefSpanToPerformAssembly, do not perform assembly
    static const unsigned MaxRefSpanToBypassAssembly = 250u;

    /// assembly parameters
    const unsigned MaxAssemblyWordSize = 76u;
    const unsigned MinAssemblyCoverage = 3u;

    /// Minimum supporting read count required to consider a haplotype for confirmation
    static const unsigned MinHaplotypeCount = 3u;

    /// If number of discovered mismatches is larger than this value,
    /// do not add mismatches into the indel buffer
    static const unsigned MaxNumMismatchesToAddToIndelBuffer = 10u;

    /// Creates an object for processing an active region
    /// \return active region object
    ActiveRegionProcessor(const known_pos_range2& posRange,
                          const pos_t prevActiveRegionEnd,
                          const reference_contig_segment& ref,
                          const unsigned maxIndelSize,
                          const unsigned sampleIndex,
                          const unsigned ploidy,
                          const GlobalAligner<int>& aligner,
                          const ActiveRegionReadBuffer& readBuffer,
                          IndelBuffer& indelBuffer,
                          CandidateSnvBuffer& candidateSnvBuffer):
        _posRange(posRange), _prevActiveRegionEnd(prevActiveRegionEnd), _ref(ref),
        _maxIndelSize(maxIndelSize), _sampleIndex(sampleIndex), _ploidy(ploidy),
        _aligner(aligner), _readBuffer(readBuffer),
        _indelBuffer(indelBuffer), _candidateSnvBuffer(candidateSnvBuffer)
    {
        _ref.get_substring(_posRange.begin_pos(), _posRange.size(), _refSegment);
    }

    /// \param pos reference position
    /// \return true if pos belongs to this active region; false otherwise
    bool contains(const pos_t pos) const
    {
        return _posRange.is_pos_intersect(pos);
    }

    /// Decompose haplotypes into primitive alleles and
    /// put them into the indel buffer and candidate SNV buffer
    void processHaplotypes();

    /// Experimental
    void addHaplotypesToExclude(const std::vector<std::string>& haplotypeToExclude);

    /// Gets selected haplotypes
    const std::vector<std::string>& getSelectedHaplotypes() const;

private:
    const known_pos_range2 _posRange;

    /// record the end position of the previous AR
    const pos_t _prevActiveRegionEnd;

    const reference_contig_segment& _ref;
    std::string _refSegment;
    const unsigned _maxIndelSize;
    const unsigned _sampleIndex;
    const unsigned _ploidy;
    const GlobalAligner<int> _aligner;

    const ActiveRegionReadBuffer& _readBuffer;

    // experimental
    std::vector<std::string> _haplotypesToExclude;

    std::vector<std::string> _selectedHaplotypes;
    std::vector<std::vector<align_id_t>> _selectedAlignIdLists;

    /// Total number of reads eligible for the haplotype generation process
    unsigned _numReadsUsedToGenerateHaplotypes;

    IndelBuffer& _indelBuffer;
    CandidateSnvBuffer& _candidateSnvBuffer;

    /// Generate haplotypes using counting and process variants
    /// \return true if haplotype generation succeeds, false otherwise
    bool generateHaplotypesWithCounting();

    /// Create haplotypes using assembly and process variants
    /// \return true if haplotype generation succeeds, false otherwise
    bool generateHaplotypesWithAssembly();

    /// Select the top haplotypes
    void selectHaplotypes(const HaplotypeToAlignIdSet& haplotypeToAlignIdSet);

    void selectOrDropHaplotypesWithSameCount(
        std::vector<std::string>& haplotypesWithSameCount,
        const HaplotypeToAlignIdSet& haplotypeToAlignIdSet,
        const bool isReferenceSelected);

    /// Do not use haplotyping to determine indel candidacy and MMDF relax positions
    void doNotUseHaplotyping();

    /// Align each selected haplotype to the reference and convert them to primitive alleles
    void processSelectedHaplotypes();

    /// Decompose haplotype into primitive alleles
    void discoverIndelsAndMismatches(
        const unsigned selectedHaplotypeIndex,
        std::vector<IndelKey>& discoveredIndelsAndMismatches,
        int& numIndels);

    /// Put discovered indels and mismatches
    /// into the indel buffer and candidate SNV buffer
    void processDiscoveredIndelsAndMismatches(
        const unsigned selectedHaplotypeIndex,
        const HaplotypeId haplotypeId,
        const std::vector<IndelKey>& discoveredIndelsAndMismatches,
        const bool doNotAddMismatchesToIndelBuffer);
};
