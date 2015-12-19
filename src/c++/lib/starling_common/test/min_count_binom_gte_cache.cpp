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
    BOOST_CHECK_EQUAL(mc.getCount(ntrials,p), min_count_binomial_gte_exact(alpha,p,ntrials));
}


BOOST_AUTO_TEST_CASE( test_mcbgc1 )
{
    testCount(1e-9,10000,1e-4);
    testCount(1e-8,100000,1e-4);
    testCount(1e-3,100000,1e-4);
}


BOOST_AUTO_TEST_SUITE_END()
