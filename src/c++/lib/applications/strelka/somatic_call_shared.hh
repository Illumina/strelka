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

/// variation on the original strawman snv caller -- implements a
/// compile-time specified grid in allele frequency space and requires
/// similar frequency as definition of non-somatic.
///

/// \author Chris Saunders
///
#ifndef __SOMATIC_CALL_SHARED_HH
#define __SOMATIC_CALL_SHARED_HH

// a simplification of diploid calls down to types relative to the reference:
//
namespace NTYPE {

enum index_t { REF,
               HOM,
               HET,
               CONFLICT,
               SIZE
             };

inline
const char*
label(const unsigned idx) {
    switch (idx) {
    case REF: return "ref";
    case HOM: return "hom";
    case HET: return "het";
    case CONFLICT: return "conflict";
    default: return "xxx";
    }
}
}

#endif
