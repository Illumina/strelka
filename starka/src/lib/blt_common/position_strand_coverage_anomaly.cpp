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
#include "blt_common/position_strand_coverage_anomaly.hh"

#include "blt_util/binomial_test.hh"



bool
position_strand_coverage_anomaly(const double alpha,
                                 const snp_pos_info& pi) {

    static const double expect_binomial_p(0.5);

    if(pi.calls.empty()) return false;

    const unsigned n_calls(pi.calls.size());

    if(n_calls<8) return false;

    unsigned n_fwd_calls(0);

    for(unsigned i(0); i<n_calls; ++i) {
        if(pi.calls[i].is_fwd_strand) n_fwd_calls += 1;
    }

    return is_reject_binomial_p(alpha,expect_binomial_p,n_fwd_calls,n_calls-n_fwd_calls);
}
