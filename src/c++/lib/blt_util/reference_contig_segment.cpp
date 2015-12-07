// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

/// \file

/// \author Chris Saunders
///

#if 0
#include "blt_util/reference_contig_segment.hh"

#include <iostream>


static
void
test_rcs()
{
    std::ostream& os(std::cerr);

    reference_contig_segment ref;
    ref.seq() = "12345";
    ref.set_offset(10);

    for (unsigned i(5); i<20; ++i)
    {
        os << "get_base i/res: " << i << " " << ref.get_base(i) << "\n";
    }

    std::string test;
    pos_t lens[2] = { 3, 8 };
    for (unsigned j(0); j<2; ++j)
    {
        for (unsigned i(5); i<20; ++i)
        {
            ref.get_substring(i,lens[j],test);
            os << "get_ss" << lens[j] << " i/res: " << i << " " << test << "\n";
        }
    }
}


//int main() { test_rcs(); }
#endif
