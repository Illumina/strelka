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

#include "rangeMap.hh"


BOOST_AUTO_TEST_SUITE( test_rangeMap )

BOOST_AUTO_TEST_CASE( test_rangeMap )
{
    rangeMap<int,int> rm;

    rm.getRef(2) = 12;
    rm.getRef(3000) = 13;
    rm.getRef(6000) = 15;

    rm.erase(2);
    rm.getRef(3000) = 3;
    rm.getRef(9000) = 12;
    rm.getRef(6000) = 2;
    rm.erase(9000);

    BOOST_REQUIRE_EQUAL(rm.getConstRef(3000), 3);
    BOOST_REQUIRE_EQUAL(rm.getConstRef(6000), 2);
}

BOOST_AUTO_TEST_CASE( test_rangeMap2 )
{
    rangeMap<int,int> rm;

    rm.getRef(10000) = 12;
    rm.getRef(9000) = 13;
    rm.getRef(7000) = 15;

    rm.erase(7000);
    BOOST_REQUIRE_EQUAL(rm.getConstRef(9000), 13);
    BOOST_REQUIRE_EQUAL(rm.getConstRef(10000), 12);
    rm.getRef(3000) = 3;
    rm.getRef(9000) = 12;
    rm.getRef(6000) = 2;
    rm.erase(9000);

    BOOST_REQUIRE_EQUAL(rm.getConstRef(3000), 3);
    BOOST_REQUIRE_EQUAL(rm.getConstRef(6000), 2);
    BOOST_REQUIRE_EQUAL(rm.getConstRef(10000), 12);
}

BOOST_AUTO_TEST_SUITE_END()

