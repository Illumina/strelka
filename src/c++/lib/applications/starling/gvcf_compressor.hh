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
 * gvcfcompressor.h
 *
 *  Created on: Feb 21, 2014
 *      Author: Morten Kallberg
 */

#pragma once

#include "gvcf_locus_info.hh"
#include "gvcf_options.hh"
#include "blt_util/RegionTracker.hh"


/// manage gVCF non-reference block compressibility criteria
struct gvcf_compressor
{
    gvcf_compressor(
        const gvcf_options& opt,
        const RegionTracker& nocompress_regions)
        : _opt(opt),
          _nocompress_regions(nocompress_regions)
    {}

    /// determine if a site should be added to current block
    /// or if the block should be written out
    bool
    is_site_compressable(
        const site_info& si) const;

    /// determine if a range of positions could be excluded
    /// (in the absence of variant signal)
    bool
    is_range_compressable(
        const known_pos_range2 range) const
    {
        if (! _opt.is_block_compression) return false;

        // check if range is in the pre-specified region that are not to be block-compressed
        return (! is_nocompress_range(range));
    }

private:

    /// is site non-compressible according to external bedfile input?:
    bool
    is_nocompress_site(
        const pos_t pos) const
    {
        return _nocompress_regions.isIntersectRegion(pos);
    }

    /// is any part of range non-compressible according to external bedfile input?:
    bool
    is_nocompress_range(
        const known_pos_range2 range) const
    {
        return _nocompress_regions.isIntersectRegion(range);
    }

    const gvcf_options& _opt;
    const RegionTracker& _nocompress_regions;
};
