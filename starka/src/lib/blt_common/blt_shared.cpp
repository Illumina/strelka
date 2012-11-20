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

/// \author Chris Saunders
///

#include "blt_common/blt_shared.hh"
#include "blt_common/position_snp_call_pprob_digt.hh"

#include "blt_util/log.hh"

#include <cstdlib>

#include <iomanip>
#include <iostream>


const char STDIN_FILENAME[] = "-";


// arbitrary... but things will be a total mess if not capped somewhere:
const unsigned MAX_FLANK_SIZE(50);



static
void
set_report_range(const blt_options& opt,
                 const pos_t ref_end,
                 pos_range& report_range){

    // The range format submitted by users on the command-line is
    // 1-indexed and fully closed, so begin=1 and end=3 refers to the
    // set: {1,2,3}
    //
    // The blt internal range format is zero-indexed and half-closed,
    // so begin=0 and end=3 refers to the set: {0,1,2} , this is the
    // internal equivilent to the command-line range specied above.
    //
    // In the next step below we change from the command-line range
    // representation to the internal representation:
    //
    report_range = opt.user_report_range;
    if(report_range.is_begin_pos){
        report_range.begin_pos -= 1;
    }

    if(! opt.is_ref_set()) return;

    if(report_range.is_begin_pos){
        if(report_range.begin_pos>=ref_end) {
            log_os << "ERROR::-report-range-begin argument must be <= reference sequence size\n";
            exit(EXIT_FAILURE);
        }
        report_range.begin_pos=std::max(report_range.begin_pos,static_cast<pos_t>(0));
    } else {
        report_range.set_begin_pos(0);
    }

    if(report_range.is_end_pos){
        report_range.end_pos=std::min(report_range.end_pos,ref_end);
    } else {
        report_range.set_end_pos(ref_end);
    }
}



static
pos_range
get_report_range_limit(const pos_range& report_range,
                       const bool is_ref_set,
                       const pos_t ref_end)
{
    pos_range rrl;

    rrl.set_begin_pos((report_range.is_begin_pos ? report_range.begin_pos : 0));
    rrl.is_end_pos = (report_range.is_end_pos || is_ref_set);
    rrl.end_pos = (report_range.is_end_pos ?
                   report_range.end_pos :
                   (is_ref_set ? ref_end : 0));
    return rrl;
}



blt_deriv_options::
blt_deriv_options(const blt_options& opt,
                  const pos_t ref_end)
    : _pdcaller(new pprob_digt_caller(opt.bsnp_diploid_theta))
{
    set_report_range(opt,ref_end,report_range);

    report_range_limit=get_report_range_limit(report_range,
                                              opt.is_ref_set(),
                                              ref_end);
}



blt_deriv_options::
~blt_deriv_options() {}



void
blt_read_counts::
report(std::ostream& os) const {
    os << "READ_COUNTS used: " << used
       << " align-score-filter: " << align_score_filter
       << " large-ref-deletion: " << large_ref_deletion
       << " unanchored: " << unanchored
       << " unmapped: " << unmapped
       << " duplicate: " << duplicate
       << " primary-analysis-filter: " << primary_filter
       << " subsample-filter: " << subsample_filter << "\n";
}
