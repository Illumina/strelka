//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
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

#include "min_count_binom_gte_cache.hh"

#include "blt_util/binomial_test.hh"



BOOST_AUTO_TEST_SUITE( mcbgc_test )


static
void
testCount(
    const double alpha,
    const unsigned ntrials,
    const double p)
{
    min_count_binom_gte_cache mc(alpha);
    BOOST_REQUIRE_EQUAL(mc.getCount(ntrials,p), min_count_binomial_gte_exact(alpha,p,ntrials));
}


BOOST_AUTO_TEST_CASE( test_count )
{
    testCount(1e-9,10000,1e-4);
    testCount(1e-8,100000,1e-4);
    testCount(1e-3,100000,1e-4);

    // test outside of Poisson range:
    testCount(1e-3,1000,0.95);
}

BOOST_AUTO_TEST_CASE( test_reject )
{
    min_count_binom_gte_cache mc(1e-9);
    BOOST_REQUIRE(! mc.isRejectNull(10000,1e-4,11));
    BOOST_REQUIRE(mc.isRejectNull(10000,1e-4,12));

    min_count_binom_gte_cache mc2(1e-3);
    BOOST_REQUIRE(! mc2.isRejectNull(100000,1e-4,11));
}

BOOST_AUTO_TEST_CASE( test_compare_reject_to_is_reject_binomial_gte_n_success_exact )
{
    static const double alpha(1e-9);
    min_count_binom_gte_cache mc(alpha);
    BOOST_REQUIRE(! mc.isRejectNull(10000,1e-4,11));
    BOOST_REQUIRE(mc.isRejectNull(10000,1e-4,12));

    BOOST_REQUIRE(! is_reject_binomial_gte_n_success_exact(alpha,1e-4,11,10000));
    BOOST_REQUIRE(is_reject_binomial_gte_n_success_exact(alpha,1e-4,12,10000));
}


BOOST_AUTO_TEST_SUITE_END()
