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

#include "somatic_indel_scoring_features.hh"

static
double
getSymmSORFeature(
    unsigned fwdAltAlleleCount,
    unsigned revAltAlleleCount,
    unsigned fwdOtherCount,
    unsigned revOtherCount)
{
    return std::log(makeSymmetric(getSampleStrandOddsRatio(fwdAltAlleleCount, revAltAlleleCount, fwdOtherCount, revOtherCount)));
}


BOOST_AUTO_TEST_SUITE( somatic_indel_scoring_features_test )


BOOST_AUTO_TEST_CASE( test_StrandOddsRatio )
{
    // test the performance of the continuous model strand bias computation at lower depths typical of a standard
    // autosome call
    static const double epsilon(0.0001);

    // no bias
    BOOST_REQUIRE_CLOSE_FRACTION(0.693147, getSymmSORFeature(0, 0, 15, 15), epsilon);

    // no bias variant
    BOOST_REQUIRE_CLOSE_FRACTION(0.693147, getSymmSORFeature(7, 7, 8, 8), epsilon);

    // all evidence on fwd strand (30x)
    BOOST_REQUIRE_CLOSE_FRACTION(0.902019, getSymmSORFeature(10, 0, 20, 0), epsilon);

    // all evidence on rev strand (30x)
    BOOST_REQUIRE_CLOSE_FRACTION(0.902019, getSymmSORFeature(0, 10, 0, 20), epsilon);

    // all evidence on fwd strand (60x)
    BOOST_REQUIRE_CLOSE_FRACTION(0.908977, getSymmSORFeature(20, 0, 40, 0), epsilon);

    // all evidence on fwd strand (30x alt 20/10)
    BOOST_REQUIRE_CLOSE_FRACTION(0.902019, getSymmSORFeature(20, 0, 10, 0), epsilon);

    // all evidence on fwd strand (30x alt 29/1)
    BOOST_REQUIRE_CLOSE_FRACTION(2.981507, getSymmSORFeature(29, 0, 1, 0), epsilon);

    // all evidence on fwd strand (30x alt 30)
    BOOST_REQUIRE_CLOSE_FRACTION(4.111142, getSymmSORFeature(30, 0, 0, 0), epsilon);

    // extreme bias
    BOOST_REQUIRE_CLOSE_FRACTION(6.867975, getSymmSORFeature(15, 0, 0, 15), epsilon);

    // moderate bias
    BOOST_REQUIRE_CLOSE_FRACTION(3.435027, getSymmSORFeature(8, 0, 8, 15), epsilon);

    // borderline bias
    BOOST_REQUIRE_CLOSE_FRACTION(2.420333, getSymmSORFeature(4, 0, 12, 15), epsilon);

    // borderline bias 2
    BOOST_REQUIRE_CLOSE_FRACTION(1.660318, getSymmSORFeature(5, 1, 10, 14), epsilon);

    // borderline bias 3
    BOOST_REQUIRE_CLOSE_FRACTION(1.157543, getSymmSORFeature(5, 2, 10, 13), epsilon);
}


BOOST_AUTO_TEST_SUITE_END()
