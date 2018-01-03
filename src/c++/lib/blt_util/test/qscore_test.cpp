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

#include "blt_util/qscore.hh"


BOOST_AUTO_TEST_SUITE( test_qscore )


BOOST_AUTO_TEST_CASE( simple_qscore_test )
{
    BOOST_REQUIRE_EQUAL(error_prob_to_qphred(0.1),10);
    BOOST_REQUIRE_EQUAL(error_prob_to_qphred(0.01),20);
    BOOST_REQUIRE_EQUAL(error_prob_to_qphred(0.001),30);
    BOOST_REQUIRE_EQUAL(error_prob_to_qphred(0.0001),40);
}


BOOST_AUTO_TEST_CASE( boundary_qscore_test )
{
    BOOST_REQUIRE_EQUAL(error_prob_to_qphred(1.0),0);
    BOOST_REQUIRE(error_prob_to_qphred(0.0) > 1000);
}


BOOST_AUTO_TEST_SUITE_END()

