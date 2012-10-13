// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
//

/// \file
///
/// \author Chris Saunders
///

#ifndef __GVCF_BLOCK_SITE_RECORD_HH
#define __GVCF_BLOCK_SITE_RECORD_HH


#include "blt_util/stream_stat.hh"
#include "starling_common/gvcf_locus_info.hh"



struct gvcf_block_site_record {

    gvcf_block_site_record(const gvcf_options& opt)
        : frac_tol(opt.block_frac_tol)
        , abs_tol(opt.block_abs_tol)
        , count(0)
    {}

    void
    reset() {
        count=0;
        block_gqx.reset();
        block_dpu.reset();
        block_dpf.reset();
    }

    // determine if the site could be joined this block:
    bool
    test(const site_info& si) const;

    // add record to this block
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



#endif
