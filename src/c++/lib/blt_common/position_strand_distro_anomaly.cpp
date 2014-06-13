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
#include "blt_common/position_strand_distro_anomaly.hh"

#include "blt_util/log.hh"
#include "blt_util/seq_util.hh"
#include "blt_util/table_test.hh"

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <iostream>


//get equiv of phred Q14.5:
const double test_max_error_prob(std::exp(1.45*std::log(10.)));



double
position_strand_distro_anomaly_pval(const snp_pos_info& pi,
#ifdef HAVE_FISHER_EXACT_TEST
                                    double* ws)
{
#else
                                    double*)
{
#endif

    static const unsigned N_STRAND(2);
    static const unsigned TABLE_SIZE(N_BASE*N_STRAND);

    const unsigned n_calls(pi.calls.size());
    int table[TABLE_SIZE];
    std::fill(table,table+TABLE_SIZE,0);
    unsigned n_good_calls(0);

    for (unsigned i(0); i<n_calls; ++i)
    {
        const uint8_t obs_id(pi.calls[i].base_id);
        assert(obs_id != BASE_ID::ANY);
        if (pi.calls[i].error_prob() > test_max_error_prob) continue;
        n_good_calls++;
        const unsigned strand(pi.calls[i].is_fwd_strand);
        table[obs_id+strand*N_BASE]++;
    }
    if (n_good_calls<2) return 1.;

#ifdef HAVE_FISHER_EXACT_TEST
    return table_exact_pval(table,N_STRAND,N_BASE,ws);
#else
    return table_chi_sqr_pval(table,N_STRAND,N_BASE);
#endif
}



bool
position_strand_distro_anomaly(const double alpha,
                               const snp_pos_info& pi,
                               double* ws)
{

    const unsigned n_calls(pi.calls.size());
    unsigned n_good_calls(0);

    for (unsigned i(0); i<n_calls; ++i)
    {
        //const char obs_base(pi.calls[i].base);
        assert(pi.calls[i].base_id != BASE_ID::ANY);
        if (pi.calls[i].error_prob() > test_max_error_prob) continue;
        n_good_calls++;
    }
    if (n_good_calls<8) return false;

    return (position_strand_distro_anomaly_pval(pi,ws) < alpha);
}
