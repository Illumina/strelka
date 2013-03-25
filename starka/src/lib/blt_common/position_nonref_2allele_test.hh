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
#ifndef __POSITION_NONREF_2ALLELE_TEST_HH
#define __POSITION_NONREF_2ALLELE_TEST_HH

#include "blt_common/blt_shared.hh"
#include "blt_common/nonref_test_call.hh"
#include "blt_common/snp_pos_info.hh"

#include <iosfwd>


namespace NR2TEST {
enum index_t {
    REF,             // reference allele only
    NONREF_MF,       // reference allele mixed with single signal allele
    NONREF_MF_NOISE, // reference allele mixed with site specific error
    NONREF_OTHER,    // reference allele mixed with non-signal allele (lhood is approximated to 0 for this state)
    SIZE
};

inline
const char*
label(const index_t i) {
    switch(i) {
    case REF: return "ref";
    case NONREF_MF: return "nonref";
    case NONREF_MF_NOISE: return "noise";
    case NONREF_OTHER: return "nonref-other";
    default: return "xxx";
    }
}
}



/// \brief Call a snp at a position under the assumption that up to
/// one nonref allele could occur at any frequency. When a snp is
/// found, optionally also provide the allele MLEs.
///
/// Method is optimized to pick up low-frequency snps, but will also
/// inevitably call higher frequency germline variation (without
/// correct priors) if present.
///
void
position_nonref_2allele_test(const snp_pos_info& pi,
                             const blt_options& opt,
                             const bool is_always_test,
                             nonref_test_call& nrc);


void
write_nonref_2allele_test(const blt_options& opt,
                          const snp_pos_info& pi,
                          const nonref_test_call& nrc,
                          std::ostream& os);

#endif
