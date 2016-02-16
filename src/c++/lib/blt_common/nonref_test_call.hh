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

/// \file

/// \author Chris Saunders
///
#ifndef __NONREF_TEST_CALL_HH
#define __NONREF_TEST_CALL_HH

#include "blt_util/seq_util.hh"

#include <boost/utility.hpp>


struct nonref_test_call : private boost::noncopyable
{

    nonref_test_call()
        : is_snp(false),
          snp_qphred(0),
          max_gt_qphred(0),
          max_gt(0),
          nonref_id(BASE_ID::ANY) {}
#if 0
    is_freq(false),
    min_loghood(0)
    {
        for (unsigned i(0); i<N_BASE; ++i) allele_freq[i]=0.;
    }
#endif

    bool is_snp;
    int snp_qphred;
    int max_gt_qphred;
    unsigned  max_gt;
    unsigned nonref_id;
#if 0
    bool is_freq;
    blt_float_t min_loghood;         // lhood from the three-parameter minimization used to find all ML frequences
    blt_float_t allele_freq[N_BASE]; // ML frequencies taking Qvalues into account.
#endif
};

#endif
