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

#include "stream_stat.hh"


BOOST_AUTO_TEST_SUITE( test_stream_stat )


BOOST_AUTO_TEST_CASE( test_stream_stat )
{
    static const double eps(0.00001);

    stream_stat ss;

    ss.add(3.);
    ss.add(4.);
    ss.add(5.);

    BOOST_REQUIRE_EQUAL(ss.size(),3);
    BOOST_REQUIRE_CLOSE(ss.min(), 3., eps);
    BOOST_REQUIRE_CLOSE(ss.max(), 5., eps);
    BOOST_REQUIRE_CLOSE(ss.mean(), 4., eps);
    BOOST_REQUIRE_CLOSE(ss.sd(), 1., eps);
}


BOOST_AUTO_TEST_SUITE_END()

