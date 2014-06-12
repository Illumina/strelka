// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
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
    switch (idx) {
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
