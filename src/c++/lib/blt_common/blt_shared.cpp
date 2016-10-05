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

///
/// \author Chris Saunders
///

#include "blt_common/blt_shared.hh"
#include "blt_common/position_snp_call_pprob_digt.hh"

#include "blt_util/log.hh"

#include <cstdlib>

#include <iomanip>
#include <iostream>
#include <sstream>


const char STDIN_FILENAME[] = "-";


// arbitrary... but things will be a total mess if not capped somewhere:
const unsigned MAX_FLANK_SIZE(100);



static
void
set_report_range(const blt_options& opt,
                 const pos_t ref_end,
                 pos_range& report_range)
{
    // The range format submitted by users on the command-line is
    // 1-indexed and fully closed, so begin=1 and end=3 refers to the
    // set: {1,2,3}
    //
    // The blt internal range format is zero-indexed and half-closed,
    // so begin=0 and end=3 refers to the set: {0,1,2} , this is the
    // internal equivalent to the command-line range specied above.
    //
    // In the next step below we change from the command-line range
    // representation to the internal representation:
    //
    report_range = opt.user_report_range;
    if (report_range.is_begin_pos)
    {
        report_range.begin_pos -= 1;
    }

    if (report_range.is_begin_pos)
    {
        if (report_range.begin_pos>=ref_end)
        {
            log_os << "ERROR::-report-range-begin argument must be <= reference sequence size\n";
            exit(EXIT_FAILURE);
        }
        report_range.begin_pos=std::max(report_range.begin_pos,static_cast<pos_t>(0));
    }
    else
    {
        report_range.set_begin_pos(0);
    }

    if (report_range.is_end_pos)
    {
        report_range.end_pos=std::min(report_range.end_pos,ref_end);
    }
    else
    {
        report_range.set_end_pos(ref_end);
    }
}



blt_deriv_options::
blt_deriv_options(
    const blt_options& opt,
    const pos_t ref_end)
    : _pdcaller(new pprob_digt_caller(opt.bsnp_diploid_theta))
{
    set_report_range(opt,ref_end,report_range);
}



/// dtor required to be in the cpp so that unique ptr can access complete data type
blt_deriv_options::
~blt_deriv_options() {}



void
blt_read_counts::
report(std::ostream& os) const
{
    os << "READ_COUNTS used: " << used
       << " align-score-filter: " << align_score_filter
       << " large-ref-deletion: " << large_ref_deletion
       << " unanchored: " << unanchored
       << " unmapped: " << unmapped
       << " duplicate: " << duplicate
       << " primary-analysis-filter: " << primary_filter
       << " subsample-filter: " << subsample_filter << "\n";
}
