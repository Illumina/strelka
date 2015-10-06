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

