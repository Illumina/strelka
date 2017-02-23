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

typedef std::vector<RangeMap<pos_t,uint8_t>> RangeSet;
typedef std::map<std::string, std::vector<align_id_t>> HaplotypeToAlignIdSet;
typedef uint8_t HaplotypeId;

/// Represent all haplotypes found in the current active region
class ActiveRegion
{
public:
    // if haplotype depth is at least HighDepth, low depth filter is not applied
    static const unsigned HighDepth = 20u;

    // if the number of reads is larger than MinNumReadsToBypassAssembly, assembly is not conducted
    static const unsigned MinNumReadsToBypassAssembly = 1000u;

    // minimum fraction of reads covering the entire region to perform counting
    float MinFracReadsCoveringRegion = 0.65f;

    // if the region length is larger than MaxRefSpanToPerformAssembly, do not perform assembly
    static const unsigned MaxRefSpanToPerformAssembly = 250u;

    // assembly parameters
    const unsigned MaxAssemblyWordSize = 76u;
    const unsigned MinAssemblyCoverage = 3u;

    // minimum haplotype count to consider
    static const unsigned MinHaplotypeCount = 3u;
    static const unsigned MaxSNVHpolSize = 4u;

    // minimum haplotype frequency to relax MMDF
    const float HaplotypeFrequencyThreshold = 0.4;

    /// Creates an active region object
    /// \param posRange position range of the active region
    /// \param ref reference
    /// \param maxIndelSize max indel size
    /// \param sampleCount sample count
    /// \param aligner aligner for aligning haplotypes to the reference
    /// \param readBuffer read buffer
    /// \return active region object
    ActiveRegion(pos_range& posRange,
                 const reference_contig_segment& ref,
                 const unsigned maxIndelSize,
                 const unsigned sampleCount,
                 const GlobalAligner<int>& aligner,
                 const ActiveRegionReadBuffer& readBuffer):
        _posRange(posRange), _ref(ref), _maxIndelSize(maxIndelSize), _sampleCount(sampleCount),
        _aligner(aligner), _readBuffer(readBuffer)
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
    bool contains(pos_t pos) const
    {
        return _posRange.is_pos_intersect(pos);
    }

    /// Decompose haplotypes into primitive alleles.
    /// Determine indel candidacy and regiser polymorphic sites to relax MMDF.
    /// \param indelBuffer
    /// \param polySites
    void processHaplotypes(IndelBuffer& indelBuffer, RangeSet& polySites);

    /// Mark a read soft-clipped
    /// \param alignId align id
    void setSoftClipped(const align_id_t alignId)
    {
        _alignIdSoftClipped.insert(alignId);
    }

    static HaplotypeId getHaplotypeId(const uint8_t rangeSetValue, const BASE_ID::index_t baseIndex)
    {
        switch (baseIndex)
        {
            case BASE_ID::A:
                return (rangeSetValue & 0x03);
            case BASE_ID::C:
                return (rangeSetValue & 0x0c) >> 2;
            case BASE_ID::G:
                return (rangeSetValue & 0x30) >> 4;
            case BASE_ID::T:
                return (rangeSetValue & 0xc0) >> 6;
            default:
                assert(false);
        }
    }

    static void addBaseId(
            const HaplotypeId complexAlleleId,
            char baseChar,
            uint8_t& rangeSetValue)
    {
        assert (complexAlleleId == 1 || complexAlleleId == 2);
        switch (baseChar)
        {
            case 'A':
                rangeSetValue |= (0x01 << (complexAlleleId-1));
                break;
            case 'C':
                rangeSetValue |= (0x04 << (complexAlleleId-1));
                break;
            case 'G':
                rangeSetValue |= (0x10 << (complexAlleleId-1));
                break;
            case 'T':
                rangeSetValue |= (0x40 << (complexAlleleId-1));
                break;
            default:
                assert(false);
        }
    }


private:
    pos_range _posRange;
    const reference_contig_segment& _ref;
    const unsigned _maxIndelSize;
    const unsigned _sampleCount;
    const GlobalAligner<int> _aligner;

    const ActiveRegionReadBuffer& _readBuffer;
    std::set<align_id_t> _alignIdSoftClipped;

    /// Create haplotypes using counting and process variants
    /// \param indelBuffer indel buffer
    /// \param polySites container to record polymorphic sites for MMDF relax
    /// \param sampleId sample id
    /// \return true if haplotype generation succeeds, false otherwise
    bool processHaplotypesWithCounting(IndelBuffer& indelBuffer, RangeSet& polySites, unsigned sampleId) const;

    /// Create haplotypes using assembly and process variants
    /// \param indelBuffer indel buffer
    /// \param polySites container to record polymorphic sites for MMDF relax
    /// \param sampleId sample id
    /// \return true if haplotype generation succeeds, false otherwise
    bool processHaplotypesWithAssembly(IndelBuffer& indelBuffer, RangeSet& polySites, unsigned sampleId) const;

    /// Bypass indels within the region
    /// \param indelBuffer indel buffer
    void bypassIndelsInBam(IndelBuffer& indelBuffer, unsigned sampleId) const;

    void convertToPrimitiveAlleles(
        const unsigned sampleId,
        const std::string& haploptypeSeq,
        const std::vector<align_id_t>& alignIdList,
        const uint8_t haplotypeId,
        IndelBuffer& indelBuffer,
        RangeSet& polySites) const;
};
