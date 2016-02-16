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

#include "boost/test/unit_test.hpp"

#include "alignment_util.hh"

//#define DEBUG_AU_TEST

#ifdef DEBUG_AU_TEST
#include <iostream>
namespace
{
std::ostream& log_os(std::cerr);
}
#endif



BOOST_AUTO_TEST_SUITE( test_alignment_util )


static
alignment
get_test_alignment()
{
    alignment al;
    al.pos=1000;
    al.is_fwd_strand=true;
    ALIGNPATH::cigar_to_apath("100M",al.path);
    return al;
}


BOOST_AUTO_TEST_CASE( test_remove_edge_deletion )
{

    alignment al(get_test_alignment());

    alignment al2(remove_edge_deletions(al));

    BOOST_CHECK_EQUAL(al, al2);

    ALIGNPATH::cigar_to_apath("3D100M",al2.path);

    alignment al3(remove_edge_deletions(al2));

    BOOST_CHECK_EQUAL(al.path, al3.path);

    BOOST_CHECK_EQUAL(static_cast<int>(al3.pos), 1003);
}



BOOST_AUTO_TEST_SUITE_END()
