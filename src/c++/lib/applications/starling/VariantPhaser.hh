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

/// \brief Short-range phasing utility using information from ActiveRegions
struct VariantPhaser : public variant_pipe_stage_base
{
    VariantPhaser(
        const starling_options& opt,
        const unsigned sampleCount,
        const std::shared_ptr<variant_pipe_stage_base>& destination)
        : variant_pipe_stage_base(destination),
          _opt(opt),
          _sampleCount(sampleCount)
    {}

    void process(std::unique_ptr<GermlineSiteLocusInfo> locusPtr) override;

    void process(std::unique_ptr<GermlineIndelLocusInfo> locusPtr) override;

private:
    const starling_options& _opt;
    const unsigned _sampleCount;

    LocusBuffer _locusBuffer;

    void flush_impl() override;

    template <class T>
    void processLocus(std::unique_ptr<T> locusPtr);

    void addLocusToBuffer(std::unique_ptr<LocusInfo> locusPtr)
    {
        _locusBuffer.push_back(std::move(locusPtr));
    }

    void outputBuffer();

    /// dump buffer contents to sink and clear object
    void createPhaseRecord(unsigned sampleIndex);

    bool
    isBuffer() const
    {
        return !(_locusBuffer.empty());
    }
};
