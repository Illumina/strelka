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

#include "fastRanksum.hh"

#include <cmath>


BOOST_AUTO_TEST_SUITE( test_fastRanksum )


BOOST_AUTO_TEST_CASE( test_fastRanksum )
{
    for (unsigned polarity(0); polarity<2; ++polarity)
    {
        const bool v1(polarity == 0);
        const bool v2(not v1);

        fastRanksum r;
        r.add_observation(v1, 44);
        r.add_observation(v1, 45); //check tie both
        r.add_observation(v2, 45);
        r.add_observation(v2, 50);
        r.add_observation(v1, 52);
        r.add_observation(v1, 53);
        r.add_observation(v1, 56);
        r.add_observation(v1, 58); //check tie single
        r.add_observation(v1, 58);
        r.add_observation(v2, 61);
        r.add_observation(v2, 63);
        r.add_observation(v1, 65);
        r.add_observation(v2, 75);
        r.add_observation(v1, 79);
        r.add_observation(v2, 85);
        r.add_observation(v2, 93);

        static const double tol(0.001);
        BOOST_REQUIRE_CLOSE(1.27021, std::abs(r.get_z_stat()), tol);
    }
}

BOOST_AUTO_TEST_CASE( test_fastRanksum2 )
{
    for (unsigned polarity(0); polarity<2; ++polarity)
    {
        const bool v1(polarity==0);
        const bool v2(not v1);
        fastRanksum r;
        r.add_observation(v1, 1);
        r.add_observation(v1, 2);
        r.add_observation(v1, 3);
        r.add_observation(v1, 4);
        r.add_observation(v1, 5);
        r.add_observation(v1, 6);
        r.add_observation(v1, 7);
        r.add_observation(v1, 8);
        r.add_observation(v1, 9);
        r.add_observation(v1, 10);
        r.add_observation(v2, 11);
        r.add_observation(v2, 12);

        static const double tol(0.001);
        BOOST_REQUIRE_CLOSE(2.14834, std::abs(r.get_z_stat()), tol);
    }
}

BOOST_AUTO_TEST_SUITE_END()

