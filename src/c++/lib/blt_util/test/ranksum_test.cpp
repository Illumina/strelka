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


BOOST_AUTO_TEST_SUITE( test_ranksum )


BOOST_AUTO_TEST_CASE( test_ranksum )
{
    ranksum r;
    r.add_observation(true,44);
    r.add_observation(true,45); //check tie both
    r.add_observation(false,45);
    r.add_observation(false,50);
    r.add_observation(true,52);
    r.add_observation(true,53);
    r.add_observation(true,56);
    r.add_observation(true,58); //check tie single
    r.add_observation(true,58);
    r.add_observation(false,61);
    r.add_observation(false,63);
    r.add_observation(true,65);
    r.add_observation(false,75);
    r.add_observation(true,79);
    r.add_observation(false,85);
    r.add_observation(false,93);

    static const double tol(0.001);
    BOOST_REQUIRE_CLOSE(-1.27021, r.get_u_stat() , tol);
}


#if 0
BOOST_AUTO_TEST_CASE( test_ranksum_runtime )
{
    ranksum r;
    for (unsigned i(0);i<10000000;++i)
    {
        r.add_observation(true,44);
        r.add_observation(true,45);
        r.add_observation(true,46);
        r.add_observation(false,45);
        r.add_observation(false,46);
        r.add_observation(false,47);
    }
    r.get_u_stat();
}
#endif


BOOST_AUTO_TEST_SUITE_END()

