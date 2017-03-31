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
struct HaplotypeIdForBase
{
    HaplotypeId a;
    HaplotypeId c;
    HaplotypeId g;
    HaplotypeId t;
};

/// for cleaning HaplotypeIdForBase
struct ClearHaplotypeIdForBase
{
    void
    operator()(HaplotypeIdForBase& val) const
    {
        val.a = 0;
        val.c = 0;
        val.t = 0;
        val.g = 0;
    }
};

typedef std::vector<RangeMap<pos_t,HaplotypeIdForBase, ClearHaplotypeIdForBase>> SampleCandidateSnvMap;

/// Stores candidate SNVs discovered in active regions
///
class CandidateSnvBuffer
{
public:
    explicit CandidateSnvBuffer(unsigned sampleCount): _candidateSnvBuffer(sampleCount)
    {}

    void addCandidateSnv(unsigned sampleId, pos_t pos, char baseChar, HaplotypeId haplotypeId);

    /// Checks if baseChar is a candidate SNV allele at this position
    /// \param sampleId sample id
    /// \param pos reference position
    /// \param baseChar read base
    /// \return true if the base matches a candidate SNV allele at this position
    bool isCandidateSnv(unsigned sampleId, pos_t pos, char baseChar) const;

    /// Gets the haplotype ID for the input single base allele
    /// \param sampleId sample id
    /// \param pos reference position
    /// \param baseIndex base index of the allele
    /// \return 0 if it's not apprearing in non-ref haplotype.
    /// 1 or 2 if it appears in one non-ref haplotype
    /// 3 if it appears in both non-ref haplotype (i.e. hetalt SNV)
    HaplotypeId getHaplotypeId(unsigned sampleId, pos_t pos, BASE_ID::index_t baseIndex) const;

    bool empty() const;

    /// Clear all candidate SNVs
    void clearSnvs();

    /// Clear all candidate SNVs up to the specified position
    /// \param pos position
    void clearUpToPos(pos_t pos);
private:
    SampleCandidateSnvMap _candidateSnvBuffer;
};