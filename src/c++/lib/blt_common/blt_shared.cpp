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

/// \file
/// \author Chris Saunders
///

#include "blt_common/blt_shared.hh"
#include "blt_common/position_snp_call_pprob_digt.hh"

#include <iostream>


// arbitrary... but things will be a total mess if not capped somewhere:
const unsigned MAX_FLANK_SIZE(100);



blt_deriv_options::
blt_deriv_options(
    const blt_options& opt)
    : _pdcaller(new pprob_digt_caller(opt.bsnp_diploid_theta))
{}



// dtor required to be in the cpp so that unique ptr can access complete data type
blt_deriv_options::
~blt_deriv_options() {}



void
blt_read_counts::
report(std::ostream& os) const
{
    os << "READ_COUNTS used: " << used
       << " align-score-filter: " << align_score_filter
       << " large-ref-deletion: " << large_ref_deletion
       << " unmapped: " << unmapped
       << " duplicate: " << duplicate
       << " primary-analysis-filter: " << primary_filter
       << " subsample-filter: " << subsample_filter << "\n";
}
