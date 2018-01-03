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

#include "algo_util.hh"


BOOST_AUTO_TEST_SUITE( test_algo_util )


BOOST_AUTO_TEST_CASE( test_getDuplicates )
{
    std::vector<int> a = {0,2,4,7,8};
    std::vector<int> b = {1,1,1,7,9};
    const auto dups(getDuplicatesInSortedInput(std::begin(a), std::end(a), std::begin(b), std::end(b)));
    BOOST_REQUIRE_EQUAL(dups.size(),1u);
    BOOST_REQUIRE_EQUAL(*(dups.begin()),7);
}


BOOST_AUTO_TEST_SUITE_END()

