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

#include "ranksum.hh"

#include <iostream>


BOOST_AUTO_TEST_SUITE( test_ranksum )


BOOST_AUTO_TEST_CASE( test_ranksum )
{
    ranksum r('A');
    r.add_observation('A',44);
    r.add_observation('A',45); //check tie both
    r.add_observation('N',45);
    r.add_observation('N',50);
    r.add_observation('A',52);
    r.add_observation('A',53);
    r.add_observation('A',56);
    r.add_observation('A',58); //check tie single
    r.add_observation('A',58);
    r.add_observation('N',61);
    r.add_observation('N',63);
    r.add_observation('A',65);
    r.add_observation('N',75);
    r.add_observation('A',79);
    r.add_observation('N',85);
    r.add_observation('N',93);

    BOOST_REQUIRE( (r.get_u_stat() < 0.00001));
}


BOOST_AUTO_TEST_SUITE_END()

