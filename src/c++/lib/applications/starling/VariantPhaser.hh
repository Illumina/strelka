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
/*
 *  Author: Sangtae Kim
 */

#pragma once

#include "gvcf_locus_info.hh"
#include "starling_shared.hh"
#include "variant_pipe_stage_base.hh"

#include "starling_common/pos_basecall_buffer.hh"

#include <iosfwd>
#include <iostream>

typedef std::vector<std::unique_ptr<LocusInfo>> LocusBuffer;

enum HaplotypeIdPair
{
    // haplotype Id 0: reference
    HomRef = 0,     // haplotype Id 0,0
    HetHap1,        // haplotype Id 0,1
    HetHap2,        // haplotype Id 0,2
    HetHap1Hap2,    // haplotype Id 1,2
    HomHap1,        // haplotype Id 1,1
    HomHap2,        // haplotype Id 2,2
    SIZE
};

/// short-range phasing utility using information from ActiveRegions
struct VariantPhaser : public variant_pipe_stage_base
{
    VariantPhaser(
        const starling_options& opt,
        const unsigned sampleCount,
        std::shared_ptr<variant_pipe_stage_base> destination)
        : variant_pipe_stage_base(destination),
          _opt(opt),
          _sampleCount(sampleCount),
          _activeRegionId(-1)
    {}

    void process(std::unique_ptr<GermlineSiteLocusInfo> locusPtr) override;

    void process(std::unique_ptr<GermlineIndelLocusInfo> locusPtr) override;

private:
    const starling_options& _opt;
    const unsigned _sampleCount;

    ActiveRegionId _activeRegionId;
    LocusBuffer _locusBuffer;

    bool _possibleHaplotypeIdPair[HaplotypeIdPair::SIZE];

    void flush_impl() override;

    template <class T>
    void processLocus(std::unique_ptr<T> locusPtr);

    void addLocusToBuffer(std::unique_ptr<LocusInfo> locusPtr)
    {
        _locusBuffer.push_back(std::move(locusPtr));
    }

    void outputBuffer();

    /// dump buffer contents to sink and clear object
    void createPhaseRecord(unsigned sampleId);

    bool
    isBuffer() const
    {
        return !(_locusBuffer.empty());
    }

    void inferHaplotypePair(const uint8_t haplotypeId, const bool isHet);
};
