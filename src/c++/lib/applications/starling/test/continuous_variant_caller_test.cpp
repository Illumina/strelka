//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
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
#include "starling_continuous_variant_caller.hh"
#include <vector>

BOOST_AUTO_TEST_SUITE( continuous_variant_caller )

BOOST_AUTO_TEST_CASE( qscore_calculation )
{
    std::vector<std::vector<unsigned>> SampleValues
    {
        //coverage,var calls, expected Q score
        {100,0,0},
        {100,1,2},
        {100,5,24},
        {200,10,43},
        {500,25,98},
        {5000,250,100}
    };

    for (const auto& item :  SampleValues)
    {
        const int FinalQValue = starling_continuous_variant_caller::getAlleleSequencingErrorQscore(item[1], item[0],
                                20, 100);
        BOOST_REQUIRE_EQUAL(FinalQValue, item[2]);
    }
}

BOOST_AUTO_TEST_CASE( compute_sb )
{
    static const double epsilon(0.0001);

    // TODO: diploid one returns -11776.1. chrM 14213
    BOOST_REQUIRE_CLOSE_FRACTION(-100.0, starling_continuous_variant_caller::strandBias(4769, 4058, 2, 4), epsilon);
    BOOST_REQUIRE_CLOSE_FRACTION(116.9888344, starling_continuous_variant_caller::strandBias(0, 62, 190, 8), epsilon);

    BOOST_REQUIRE_CLOSE_FRACTION(251.724361, starling_continuous_variant_caller::strandBias(379, 0, 4940, 5601), epsilon);
    BOOST_REQUIRE_CLOSE_FRACTION(-100.0, starling_continuous_variant_caller::strandBias(379, 400, 4940, 5601), epsilon);
}


BOOST_AUTO_TEST_SUITE_END()
