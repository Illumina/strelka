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
///
/// \author Chris Saunders
///

#ifndef __REPORT_BACON_CALLS_HH
#define __REPORT_BACON_CALLS_HH

#include "blt_common/blt_shared.hh"
#include "blt_common/blt_streams.hh"
#include "blt_common/position_snp_call_bacon.hh"

#include <boost/utility.hpp>



enum bacon_snp_t {
    BACON_SNP_DIFF,
    BACON_SNP_HET1,
    BACON_SNP_HET2,
    BACON_SNP_HETO
};



struct bacon_info : private boost::noncopyable {

    bacon_info()
        : unused_count(0), used_count(0),
          is_snp(false), is_het(false), is_max2_reported(false) {
        for(unsigned i(0);i<N_BASE;++i) base_count[i] = 0;
    }

    bacon_scores bas;
    unsigned unused_count;
    unsigned used_count;
    unsigned base_count[N_BASE];
    bool is_snp;
    bool is_het;
    bool is_max2_reported;
    bacon_snp_t bst;
};


void
get_bacon_scores(const blt_options& client_opt,
                 const snp_pos_info& pi,
                 const unsigned unused_count,
                 bacon_info& bi);

void
report_bacon_allele_call(const blt_streams& client_io,
                         const unsigned pos,
                         const bacon_info& bi);

void
report_bacon_snp_call(const blt_streams& client_io,
                      const unsigned pos,
                      const char ref_base,
                      const bacon_info& bi);

#endif
