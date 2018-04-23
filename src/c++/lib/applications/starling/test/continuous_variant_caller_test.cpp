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

#include "starling_continuous_variant_caller.hh"

#include <vector>

BOOST_AUTO_TEST_SUITE( continuous_variant_caller_test_suite )

BOOST_AUTO_TEST_CASE( qscore_calculation_test )
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

BOOST_AUTO_TEST_CASE( computeContinuousModelStrandBiasTest )
{
    // several strand bias tests exercising cases from mitochondria
    static const double epsilon(0.0001);

    // TODO: diploid one returns -11776.1. chrM 14213
    BOOST_REQUIRE(starling_continuous_variant_caller::strandBias(4769, 4058, 2, 4) < -100);
    BOOST_REQUIRE_CLOSE_FRACTION(116.9888344, starling_continuous_variant_caller::strandBias(0, 62, 190, 8), epsilon);

    BOOST_REQUIRE_CLOSE_FRACTION(251.724361, starling_continuous_variant_caller::strandBias(379, 0, 4940, 5601), epsilon);
    BOOST_REQUIRE(starling_continuous_variant_caller::strandBias(379, 400, 4940, 5601) < -100);
}


BOOST_AUTO_TEST_CASE( computeContinuousModelStrandBiasLowDepthTest )
{
    // test the performance of the continuous model strand bias computation at lower depths typical of a standard
    // autosome call
    static const double epsilon(0.0001);

    // no bias
    BOOST_REQUIRE_CLOSE_FRACTION(-0.075188, starling_continuous_variant_caller::strandBias(0, 0, 15, 15), epsilon);

    // no bias variant
    BOOST_REQUIRE_CLOSE_FRACTION(-26.764472, starling_continuous_variant_caller::strandBias(7, 7, 8, 8), epsilon);

    // all evidence on fwd strand (30x)
    BOOST_REQUIRE_CLOSE_FRACTION(0, starling_continuous_variant_caller::strandBias(10, 0, 20, 0), epsilon);

    // all evidence on rev strand (30x)
    BOOST_REQUIRE_CLOSE_FRACTION(0, starling_continuous_variant_caller::strandBias(0, 10, 0, 20), epsilon);

    // all evidence on fwd strand (60x)
    BOOST_REQUIRE_CLOSE_FRACTION(0, starling_continuous_variant_caller::strandBias(20, 0, 40, 0), epsilon);

    // all evidence on fwd strand (30x alt 20/10)
    BOOST_REQUIRE_CLOSE_FRACTION(0, starling_continuous_variant_caller::strandBias(20, 0, 10, 0), epsilon);

    // all evidence on fwd strand (30x alt 29/1)
    BOOST_REQUIRE_CLOSE_FRACTION(0, starling_continuous_variant_caller::strandBias(29, 0, 1, 0), epsilon);

    // all evidence on fwd strand (30x alt 30)
    BOOST_REQUIRE_CLOSE_FRACTION(0, starling_continuous_variant_caller::strandBias(30, 0, 0, 0), epsilon);

    // extreme bias
    BOOST_REQUIRE_CLOSE_FRACTION(20.719227, starling_continuous_variant_caller::strandBias(15, 0, 0, 15), epsilon);

    // moderate bias
    BOOST_REQUIRE_CLOSE_FRACTION(6.536161, starling_continuous_variant_caller::strandBias(8, 0, 8, 15), epsilon);

    // borderline bias
    BOOST_REQUIRE_CLOSE_FRACTION(2.848280, starling_continuous_variant_caller::strandBias(4, 0, 12, 15), epsilon);

    // borderline bias 2
    BOOST_REQUIRE_CLOSE_FRACTION(0.095867, starling_continuous_variant_caller::strandBias(5, 1, 10, 14), epsilon);

    // borderline bias 2 (60x)
    BOOST_REQUIRE_CLOSE_FRACTION(0.191734, starling_continuous_variant_caller::strandBias(10, 2, 20, 28), epsilon);

    // borderline bias 3
    BOOST_REQUIRE_CLOSE_FRACTION(-3.911326, starling_continuous_variant_caller::strandBias(5, 2, 10, 13), epsilon);

    // borderline bias 3 (60x)
    BOOST_REQUIRE_CLOSE_FRACTION(-7.822653, starling_continuous_variant_caller::strandBias(10, 4, 20, 26), epsilon);
}

BOOST_AUTO_TEST_SUITE_END()
