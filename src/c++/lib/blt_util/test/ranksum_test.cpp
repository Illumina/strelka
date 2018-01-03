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
    BOOST_REQUIRE_CLOSE(-1.27021, r.get_u_stat(), tol);
}


#if 0
BOOST_AUTO_TEST_CASE( test_ranksum_runtime )
{
    ranksum r;
    for (unsigned i(0); i<10000000; ++i)
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

