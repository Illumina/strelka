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
 *  Author: Morten Kallberg
 */

#include "gvcf_compressor.hh"

#include "blt_util/RegionTracker.hh"



gvcf_compressor::
gvcf_compressor(
    const gvcf_options& opt,
    const RegionTracker& nocompress_regions)
    : _opt(opt),
      _nocompress_regions(nocompress_regions)
{}



bool
gvcf_compressor::
is_site_compressable(
    const GermlineSiteCallInfo& si) const
{
    if (si.forcedOutput) return false;
    if (! _opt.is_block_compression) return false;

    if (si.is_snp()) return false;

    if ((si.ref != 'N') && (si.n_used_calls > 0))
    {
        const double reffrac(static_cast<double>(si.alleleObservationCounts(base_to_id(si.ref))) /
                             static_cast<double>(si.n_used_calls));
        if ((reffrac+_opt.block_max_nonref) <= 1) return false;
    }

    // check if site is in the pre-specified region that are not to be block-compressed
    return (! _nocompress_regions.isIntersectRegion(si.pos));
}



bool
gvcf_compressor::
is_range_compressable(
    const known_pos_range2 range) const
{
    if (! _opt.is_block_compression) return false;

    // check if range is in the pre-specified region that are not to be block-compressed
    return (! _nocompress_regions.isIntersectRegion(range));
}
