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
    switch(idx) {
    case REF: return "ref";
    case HOM: return "hom";
    case HET: return "het";
    case CONFLICT: return "conflict";
    default: return "xxx";
    }
}
}

#endif
