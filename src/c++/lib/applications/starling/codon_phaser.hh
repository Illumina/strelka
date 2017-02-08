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
 *  Test for codon-phasing.
 *
 *  Created on: Aug 10, 2013
 *  Author: Morten Kallberg
 */

#pragma once

#include "gvcf_locus_info.hh"
#include "starling_shared.hh"
#include "variant_pipe_stage_base.hh"

#include "starling_common/pos_basecall_buffer.hh"

#include <iosfwd>

typedef std::vector<std::unique_ptr<GermlineSiteLocusInfo>> SiteLocusBuffer;
typedef std::vector<std::unique_ptr<GermlineIndelLocusInfo>> IndelLocusBuffer;

/// short-range phasing utility for het-snps
///
/// requires extended preservation of the pileup buffer so that it
/// can go back and recover phase information form a candidate phasing block
///
/// \TODO generally recognized development direction is to record some kind of
///       read id in SNP pileups and indel support info so that we can go back
///       and phase from the hets without having to keep the whole read buffer (and so
///       read filtration, etc. is an exact match to the pileup).
///       Will this be worth doing before we transition to a haplotype assembly model
///       for short-range phasing?
///
struct Codon_phaser : public variant_pipe_stage_base
{
    Codon_phaser(
        const starling_options& opt,
        const unsigned sampleCount,
        std::shared_ptr<variant_pipe_stage_base> destination)
        : variant_pipe_stage_base(destination),
          _opt(opt),
          _sampleCount(sampleCount),
          _activeRegionId(-1),
          _sampleSiteLocusBuffer(sampleCount),
          _sampleIndelLocusBuffer(sampleCount)
    {}

    void process(std::unique_ptr<GermlineSiteLocusInfo> locusPtr) override;

    void process(std::unique_ptr<GermlineIndelLocusInfo> locusPtr) override;

    void create_phased_record(const unsigned sampleIndex);

    bool
    isBuffer(unsigned sampleId) const
    {
        return !(_sampleSiteLocusBuffer[sampleId].empty() and _sampleIndelLocusBuffer[sampleId].empty());
    }

private:
    void flush_impl() override;

    static
    bool
    isPhasableLocus(
            const LocusInfo &locus,
            const unsigned sampleIndex)
    {
        if (not locus.isVariantLocus()) return false;
        return locus.getSample(sampleIndex).max_gt().isHet();
    }

    void addSiteLocusToBuffer(std::unique_ptr<GermlineSiteLocusInfo> locusPtr)
    {
        for (unsigned sampleId(0); sampleId < _sampleCount; ++sampleId)
        {
            if (isPhasableLocus(*locusPtr, sampleId))
            {
                _sampleSiteLocusBuffer[sampleId].push_back(std::move(locusPtr));
            }
        }
    }

    void addIndelLocusToBuffer(std::unique_ptr<GermlineIndelLocusInfo> locusPtr)
    {
        for (unsigned sampleId(0); sampleId < _sampleCount; ++sampleId)
        {
            if (isPhasableLocus(*locusPtr, sampleId))
            {
                _sampleIndelLocusBuffer[sampleId].push_back(std::move(locusPtr));
            }
        }
    }

    void outputBuffer();

    /// dump buffer contents to sink and clear object
    void outputBuffer(unsigned sampleId);

    const starling_options& _opt;
    const unsigned _sampleCount;

    ActiveRegionId _activeRegionId;
    std::vector<SiteLocusBuffer> _sampleSiteLocusBuffer;
    std::vector<IndelLocusBuffer> _sampleIndelLocusBuffer;
};
