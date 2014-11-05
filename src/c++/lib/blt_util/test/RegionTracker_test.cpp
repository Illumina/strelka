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

#define private public

#include "RegionTracker.hh"


BOOST_AUTO_TEST_SUITE( test_RegionTracker )


BOOST_AUTO_TEST_CASE( test_RegionTracker )
{
    // simplist test
    RegionTracker rt;

    rt.addRegion(known_pos_range2(0,1));
    BOOST_REQUIRE(rt.isInRegion(0));
    BOOST_REQUIRE(! rt.isInRegion(1));
}


BOOST_AUTO_TEST_CASE( test_RegionTracker2 )
{
    // region overlap tests
    {
        RegionTracker rt;

        rt.addRegion(known_pos_range2(5,10));
        rt.addRegion(known_pos_range2(2,3));
        BOOST_REQUIRE_EQUAL(rt._regions.size(),2);
        BOOST_REQUIRE(  rt.isInRegion(2));
        BOOST_REQUIRE(! rt.isInRegion(3));
        BOOST_REQUIRE(! rt.isInRegion(4));
        BOOST_REQUIRE(  rt.isInRegion(5));

        rt.addRegion(known_pos_range2(3,7));
        BOOST_REQUIRE_EQUAL(rt._regions.size(),1);
        BOOST_REQUIRE(  rt.isInRegion(4));
    }
    {
        RegionTracker rt;
        rt.addRegion(known_pos_range2(5,10));
        rt.addRegion(known_pos_range2(2,3));
        rt.addRegion(known_pos_range2(2,5));
        BOOST_REQUIRE_EQUAL(rt._regions.size(),1);
        BOOST_REQUIRE(  rt.isInRegion(4));
    }

    {
        RegionTracker rt;
        rt.addRegion(known_pos_range2(5,10));
        rt.addRegion(known_pos_range2(2,3));
        rt.addRegion(known_pos_range2(2,4));
        BOOST_REQUIRE_EQUAL(rt._regions.size(),2);
        BOOST_REQUIRE(  rt.isInRegion(3));
    }

    {
        RegionTracker rt;
        rt.addRegion(known_pos_range2(5,10));
        rt.addRegion(known_pos_range2(2,3));
        rt.addRegion(known_pos_range2(4,5));
        BOOST_REQUIRE_EQUAL(rt._regions.size(),2);
        BOOST_REQUIRE(  rt.isInRegion(4));
    }

    {
        RegionTracker rt;
        rt.addRegion(known_pos_range2(1,10));
        rt.addRegion(known_pos_range2(4,5));
        BOOST_REQUIRE_EQUAL(rt._regions.size(),1);
        BOOST_REQUIRE(  rt.isInRegion(4));
    }

    {
        RegionTracker rt;
        rt.addRegion(known_pos_range2(4,5));
        rt.addRegion(known_pos_range2(1,10));
        BOOST_REQUIRE_EQUAL(rt._regions.size(),1);
        BOOST_REQUIRE(  rt.isInRegion(4));
    }
}


BOOST_AUTO_TEST_CASE( test_RegionTracker3 )
{
    // region remove tests
    RegionTracker rt;

    rt.addRegion(known_pos_range2(5,10));
    rt.addRegion(known_pos_range2(2,3));
    rt.addRegion(known_pos_range2(14,15));
    rt.addRegion(known_pos_range2(24,25));
    BOOST_REQUIRE_EQUAL(rt._regions.size(),4);

    rt.removeToPos(2);
    BOOST_REQUIRE_EQUAL(rt._regions.size(),3);
    rt.removeToPos(2);
    BOOST_REQUIRE_EQUAL(rt._regions.size(),3);
    rt.removeToPos(6);
    BOOST_REQUIRE_EQUAL(rt._regions.size(),3);
    rt.removeToPos(16);
    BOOST_REQUIRE_EQUAL(rt._regions.size(),1);
}


BOOST_AUTO_TEST_SUITE_END()

