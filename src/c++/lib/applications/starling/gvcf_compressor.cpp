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
 * gvcfcompressor.cpp
 *
 *  Created on: Feb 21, 2014
 *  Author: Morten Kallberg
 */

#include "gvcf_compressor.hh"



#if 0
/// determine the maximum no-call range that can be compressed, return the number of loci the block can be extended by
///
/// appears to be unused? -- note if this is needed, a fast version could be built into RegionTracker
int gvcf_compressor::max_compressible_nocall_range(const int start, const int end)
{
    if (!this->minor_allele_loaded)
        return (end-start);
    int i;
    for (i=start; i<=end; i++)
        if (this->is_nocompress_site(i))
            return (i-start);

    return (end-start);
}
#endif



bool
gvcf_compressor::
is_site_compressable(
    const site_info& si) const
{
    if (si.forcedOutput) return false;
    if (! _opt.is_block_compression) return false;

    if (si.is_snp()) return false;

    if (si.ref!='N' && si.n_used_calls > 0)
    {
        const double reffrac(static_cast<double>(si.known_counts[base_to_id(si.ref)]) /
                             static_cast<double>(si.n_used_calls));
        if (reffrac+_opt.block_max_nonref <= 1) return false;
    }

    // check if site is in the pre-specified region that are not to be block-compressed
    return (! is_nocompress_site(si.pos));
}
