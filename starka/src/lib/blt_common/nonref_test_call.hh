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
#ifndef __NONREF_TEST_CALL_HH
#define __NONREF_TEST_CALL_HH

#include "blt_util/seq_util.hh"

#include <boost/utility.hpp>


struct nonref_test_call : private boost::noncopyable {

    nonref_test_call()
        : is_snp(false),
          snp_qphred(0),
          max_gt_qphred(0),
          max_gt(0),
          nonref_id(BASE_ID::ANY) {}
#if 0
    is_freq(false),
    min_loghood(0) {
        for(unsigned i(0); i<N_BASE; ++i) allele_freq[i]=0.;
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
