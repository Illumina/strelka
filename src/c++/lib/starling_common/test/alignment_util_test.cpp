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
