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

///
/// \author Chris Saunders
///

#pragma once


#include "gvcf_locus_info.hh"
#include "gvcf_options.hh"
#include "blt_util/stream_stat.hh"


/// manages compressed site record blocks output in the gVCF
///
struct gvcf_block_site_record
{
    gvcf_block_site_record(const gvcf_options& opt)
        : frac_tol(static_cast<double>(opt.block_percent_tol)/100.)
        , abs_tol(opt.block_abs_tol)
        , count(0)
    {}

    void
    reset()
    {
        count=0;
        block_gqx.reset();
        block_dpu.reset();
        block_dpf.reset();
        record.smod.clear(); //clear filter as well
    }

    /// determine if the given site could be joined to this block:
    bool
    test(const site_info& si) const;

    /// add site to the current block
    void
    join(const site_info& si);

    site_info record;

    const double frac_tol;
    const int abs_tol;
    int count;
    stream_stat block_gqx;
    stream_stat block_dpu;
    stream_stat block_dpf;
    //stream_stat _blockMQ;
};
