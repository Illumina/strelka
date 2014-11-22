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

#include "blt_util/binomial_test.hh"



BOOST_AUTO_TEST_SUITE( test_binomial_test )


BOOST_AUTO_TEST_CASE( test_simple_binomial_test )
{
    static const double alpha(0.01);
    BOOST_REQUIRE(! is_reject_binomial_p(alpha,0.5,1,10));
    BOOST_REQUIRE(  is_reject_binomial_p(alpha,0.5,10,100));

    // run the counts high enough to hit the chi-sq switchpoint:
    BOOST_REQUIRE(  is_reject_binomial_p(alpha,0.5,100,1000));
    BOOST_REQUIRE(! is_reject_binomial_p(alpha,0.1,100,1000));
}


BOOST_AUTO_TEST_SUITE_END()
