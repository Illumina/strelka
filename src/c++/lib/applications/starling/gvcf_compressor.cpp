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
    if (! _opt.is_block_compression) return false;

    if (si.dgt.is_snp) return false;

    if (si.ref!='N')
    {
        const double reffrac(static_cast<double>(si.known_counts[si.dgt.ref_gt]) /
                             static_cast<double>(si.n_used_calls));
        if (reffrac+_opt.block_max_nonref <= 1) return false;
    }

    // check if site is in the pre-specified region that are not to be block-compressed
    return (! is_nocompress_site(si.pos));
}
