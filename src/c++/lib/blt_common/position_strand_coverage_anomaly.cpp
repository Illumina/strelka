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

/// \file

/// \author Chris Saunders
///
#include "blt_common/position_strand_coverage_anomaly.hh"

#include "blt_util/binomial_test.hh"



bool
position_strand_coverage_anomaly(const double alpha,
                                 const snp_pos_info& pi) {

    static const double expect_binomial_p(0.5);

    if (pi.calls.empty()) return false;

    const unsigned n_calls(pi.calls.size());

    if (n_calls<8) return false;

    unsigned n_fwd_calls(0);

    for (unsigned i(0); i<n_calls; ++i) {
        if (pi.calls[i].is_fwd_strand) n_fwd_calls += 1;
    }

    return is_reject_binomial_p(alpha,expect_binomial_p,n_fwd_calls,n_calls-n_fwd_calls);
}
