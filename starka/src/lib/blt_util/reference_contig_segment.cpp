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

#include "blt_util/reference_contig_segment.hh"

#include <iostream>



void
test_rcs() {
    std::ostream& os(std::cerr);

    reference_contig_segment ref;
    ref.seq() = "12345";
    ref.set_offset(10);

    for(unsigned i(5); i<20; ++i) {
        os << "get_base i/res: " << i << " " << ref.get_base(i) << "\n";
    }

    std::string test;
    pos_t lens[2] = { 3, 8 };
    for(unsigned j(0); j<2; ++j) {
        for(unsigned i(5); i<20; ++i) {
            ref.get_substring(i,lens[j],test);
            os << "get_ss" << lens[j] << " i/res: " << i << " " << test << "\n";
        }
    }
}


//int main() { test_rcs(); }

