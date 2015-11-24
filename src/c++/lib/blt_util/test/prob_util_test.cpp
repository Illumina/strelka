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

#include "prob_util.hh"

BOOST_AUTO_TEST_SUITE( prob_util )

BOOST_AUTO_TEST_CASE( test_prob_comp )
{
    static const float eps = 0.00000001;

    // prob simulates a typical strong call posterior where we have strong
    // evidence for state 0 and expect a high q-score:
    const float prob[] = {0.99999999,0.000000002,0.000000003,0.000000005};
    const unsigned psize(sizeof(prob)/sizeof(float));

    const float expect = 0.00000001;
    BOOST_REQUIRE_CLOSE(prob_comp(prob,prob+psize,0), expect, eps);

    // uncomment this case to demo why prob_comp is used:
#if 0
    const float val1(1.-prob[0]);
    BOOST_REQUIRE_CLOSE(val1, expect , eps);
#endif
}

BOOST_AUTO_TEST_SUITE_END()

