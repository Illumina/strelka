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

#include "blt_util/reference_contig_segment.hh"

#include <iostream>



static
void
test_rcs() {
    std::ostream& os(std::cerr);

    reference_contig_segment ref;
    ref.seq() = "12345";
    ref.set_offset(10);

    for(unsigned i(5);i<20;++i) {
        os << "get_base i/res: " << i << " " << ref.get_base(i) << "\n";
    }

    std::string test;
    pos_t lens[2] = { 3, 8 };
    for(unsigned j(0);j<2;++j) {
        for(unsigned i(5);i<20;++i) {
            ref.get_substring(i,lens[j],test);
            os << "get_ss" << lens[j] << " i/res: " << i << " " << test << "\n";
        }
    }
}


//int main() { test_rcs(); }

