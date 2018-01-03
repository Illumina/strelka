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
 *      Author: Morten Kallberg
 */

#pragma once

#include "gvcf_locus_info.hh"
#include "gvcf_options.hh"

// fwd this class b/c it's an iostream polluter
struct RegionTracker;


/// manage gVCF non-reference block compression criteria
struct gvcf_compressor
{
    gvcf_compressor(
        const gvcf_options& opt,
        const RegionTracker& nocompress_regions);

    /// determine if a site should be added to current block
    /// or if the block should be written out
    bool
    is_site_compressible(
        const GermlineSiteLocusInfo& locus) const;

    /// determine if a range of positions could be excluded
    /// (in the absence of variant signal)
    bool
    is_range_compressible(
        const known_pos_range2 range) const;

private:
    const gvcf_options& _opt;
    const RegionTracker& _nocompress_regions;
};
