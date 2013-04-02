// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

/// \file

/// \author Chris Saunders
///
#ifndef __POSITION_NONREF_TEST_HH
#define __POSITION_NONREF_TEST_HH

#include "blt_common/blt_shared.hh"
#include "blt_common/nonref_test_call.hh"
#include "blt_common/snp_pos_info.hh"

#include <iosfwd>


namespace NRTEST {
enum index_t {
    REF,
    NONREF,
    SIZE
};

inline
const char*
label(const index_t i) {
    switch(i) {
    case REF: return "ref";
    case NONREF: return "nonref";
    default: return "xxx";
    }
}
}


/// \brief Call a snp at a position under the assumption that any
/// combination of non-reference requencies could occur. When a snp is
/// found, optionally also provide the allele MLEs.
///
void
position_nonref_test(const snp_pos_info& pi,
                     const double nonref_site_prob,
                     const double min_nonref_freq,
                     const bool is_mle_freq,
                     nonref_test_call& nrc);


void
write_nonref_test(const blt_options& opt,
                  const snp_pos_info& pi,
                  const nonref_test_call& nrc,
                  std::ostream& os);

#endif
