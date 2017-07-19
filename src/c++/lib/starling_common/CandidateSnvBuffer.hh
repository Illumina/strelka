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
    explicit CandidateSnvBuffer()
    {}

    void addCandidateSnv(
        const pos_t pos,
        const char baseChar,
        const HaplotypeId haplotypeId,
        const float altHaplotypeCountRatio);

    /// Checks if baseChar is a candidate SNV allele at this position
    /// \param pos reference position
    /// \param baseChar read base
    /// \return true if the base matches a candidate SNV allele at this position
    bool isCandidateSnv(const pos_t pos, const char baseChar) const;

    /// Gets the haplotype ID for the input single base allele
    /// \param pos reference position
    /// \param baseIndex base index of the allele
    /// \return 0 if it's not appearing in non-ref haplotype.
    /// 1 or 2 if it appears in one non-ref haplotype
    /// 3 if it appears in both non-ref haplotype (i.e. hetalt SNV)
    HaplotypeId getHaplotypeId(const pos_t pos, const BASE_ID::index_t baseIndex) const;

    float getAltHaplotypeCountRatio(const pos_t pos) const;

    bool empty() const;

    /// Clear all candidate SNVs up to the specified position
    /// \param pos position
    void clearUpToPos(const pos_t pos);
private:
    SampleCandidateSnvMap _candidateSnvBuffer;
};
