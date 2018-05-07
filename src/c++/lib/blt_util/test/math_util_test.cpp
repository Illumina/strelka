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

#include "math_util.hh"


BOOST_AUTO_TEST_SUITE( math_util_test_suite )

static
void
medianChecker(
    const int trueVal,
    const std::vector<int>& nums)
{
    std::vector<int> copy1(nums);
    BOOST_REQUIRE_EQUAL(trueVal,_ps_median(std::begin(copy1), std::end(copy1)));
    std::vector<int> copy2(nums);
    BOOST_REQUIRE_EQUAL(trueVal,_ne_median(std::begin(copy2), std::end(copy2)));
}


BOOST_AUTO_TEST_CASE( test_median )
{
    {
        std::vector<int> nums = {9,8,7,0,2,1,3};
        medianChecker(3,nums);
    }
    {
        std::vector<int> nums  = {9,8,7,0,2,1,3,10,11};
        medianChecker(7,nums);
    }
}


BOOST_AUTO_TEST_SUITE_END()

