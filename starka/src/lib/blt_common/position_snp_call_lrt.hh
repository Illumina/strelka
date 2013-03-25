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
#ifndef __POSITION_SNP_CALL_LRT_HH
#define __POSITION_SNP_CALL_LRT_HH

#include "blt_common/blt_shared.hh"

#include <boost/utility.hpp>

#include <iosfwd>


struct lrt_snp_call : private boost::noncopyable {

    lrt_snp_call()
        : is_snp(false),
          null_loghood(0),
          min_test_loghood(0),
          min_loghood(0),
          snp_prob(0) {
        for(unsigned i(0); i<N_BASE; ++i) allele_freq[i]=0.;
    }

    bool is_snp;
    double null_loghood;
    double min_test_loghood; // lhood from the single-parameter test used to determine sigificance
    double min_loghood;      // lhood from the three-parameter minimization used to find all ML frequences
    double snp_prob;
    double allele_freq[N_BASE]; // ML frequencies taking Qvalues into account.
};

std::ostream& operator<<(std::ostream& os,
                         const lrt_snp_call& sc);



/// \brief call a snp @ pos by likelihood ratio test
///
/// just returns yes/no for a snp right now, this will have to return
/// some more useful info eventually...
///
void
position_snp_call_lrt(const double alpha,
                      const snp_pos_info& pi,
                      lrt_snp_call& sc);

#endif
