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

///
/// \author Chris Saunders
///
#include "blt_common/position_strand_coverage_anomaly.hh"

#include "blt_util/binomial_test.hh"



bool
position_strand_coverage_anomaly(const double alpha,
                                 const snp_pos_info& pi)
{
    static const double expect_binomial_p(0.5);

    if (pi.calls.empty()) return false;

    const unsigned n_calls(pi.calls.size());

    if (n_calls<8) return false;

    unsigned n_fwd_calls(0);

    for (unsigned i(0); i<n_calls; ++i)
    {
        if (pi.calls[i].is_fwd_strand) n_fwd_calls += 1;
    }

    return is_reject_binomial_twosided(alpha,expect_binomial_p,n_fwd_calls,n_calls);
}
