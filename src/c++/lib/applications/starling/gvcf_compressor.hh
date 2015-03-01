// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
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
