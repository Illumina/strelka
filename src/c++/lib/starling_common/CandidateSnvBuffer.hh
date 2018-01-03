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

#include <cstdint>
#include <blt_util/blt_types.hh>
#include <blt_util/seq_util.hh>
#include <blt_util/RangeMap.hh>

typedef uint8_t HaplotypeId;

/// stores haplotype id for A, C, T, G
struct HaplotypeIdAndCountRatio
{
    HaplotypeId a;
    HaplotypeId c;
    HaplotypeId g;
    HaplotypeId t;
    float altHaplotypeCountRatio;
};

/// for cleaning HaplotypeIdForBase
struct ClearHaplotypeIdForBase
{
    void
    operator()(HaplotypeIdAndCountRatio& val) const
    {
        val.a = 0;
        val.c = 0;
        val.t = 0;
        val.g = 0;
        val.altHaplotypeCountRatio = 0.0f;
    }
};

typedef RangeMap<pos_t,HaplotypeIdAndCountRatio, ClearHaplotypeIdForBase> SampleCandidateSnvMap;

/// Stores candidate SNVs discovered in active regions
///
class CandidateSnvBuffer
{
public:
    explicit CandidateSnvBuffer(const unsigned sampleCount):
        _sampleCount(sampleCount),
        _candidateSnvBuffer(sampleCount)
    {}

    void addCandidateSnv(
        const unsigned sampleIndex,
        const pos_t pos,
        const char baseChar,
        const HaplotypeId haplotypeId,
        const float altHaplotypeCountRatio);

    /// Checks if baseChar is a candidate SNV allele at this position in any sample
    /// \return true if the base matches a candidate SNV allele at this position in any sample
    bool isCandidateSnvAnySample(const pos_t pos, const char baseChar) const;

    /// Gets the haplotype ID for the input single base allele
    /// \param sampleIndex sample index
    /// \param pos reference position
    /// \param baseIndex base index of the allele
    /// \return 0 if it's not appearing in non-ref haplotype.
    /// 1 or 2 if it appears in one non-ref haplotype
    /// 3 if it appears in both non-ref haplotype (i.e. hetalt SNV)
    HaplotypeId getHaplotypeId(const unsigned sampleIndex, const pos_t pos, const BASE_ID::index_t baseIndex) const;

    float getAltHaplotypeCountRatio(const unsigned sampleIndex, const pos_t pos) const;

    /// Returns true if no candidate SNV exists in any sample
    bool empty() const;

    /// Clear all candidate SNVs up to the specified position
    /// \param sampleIndex sample index
    /// \param pos position
    void clearUpToPos(const unsigned sampleIndex, const pos_t pos);

    void clearSnvs()
    {
        for (unsigned sampleIndex(0); sampleIndex<_sampleCount; ++sampleIndex)
            _candidateSnvBuffer[sampleIndex].clear();
    }
private:
    unsigned _sampleCount;
    std::vector<SampleCandidateSnvMap> _candidateSnvBuffer;
};
