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
#ifndef __MONOGT_HH
#define __MONOGT_HH


namespace MONOGT {
enum index_t {
    A,
    C,
    G,
    T,
    SIZE
};

inline
const char*
label(const unsigned idx) {
    switch(idx) {
    case A: return "A";
    case C: return "C";
    case G: return "G";
    case T: return "T";
    default: return "X";
    }
}

// the lhood function is no longer so general that these values can actually be changed...
inline
double
expect(const int base_id,
       const int gt) {

    static const unsigned N_BASE(4);

    static const double ex[SIZE][N_BASE] = {{ 1.0, 0.0, 0.0, 0.0},
        { 0.0, 1.0, 0.0, 0.0},
        { 0.0, 0.0, 1.0, 0.0},
        { 0.0, 0.0, 0.0, 1.0}
    };

    return ex[gt][base_id];
}
}


#endif
