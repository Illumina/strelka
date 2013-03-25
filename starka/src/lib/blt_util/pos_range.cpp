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

#include "blt_util/pos_range.hh"

#include <iostream>


// output is always 1-indexed inclusive interval:
//
std::ostream& operator<<(std::ostream& os, const pos_range& pr) {

    os << "[";
    if(pr.is_begin_pos) {
        os << pr.begin_pos+1;
    } else {
        os << "-inf";
    }
    os << " .. ";
    if(pr.is_end_pos) {
        os << pr.end_pos;
    } else {
        os << "inf";
    }
    os << "]";

    return os;
}
