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

///
/// \author Eunho Noh
///

#include "calibration/IndelErrorModel.hh"
#include "boost/test/unit_test.hpp"

BOOST_AUTO_TEST_SUITE( test_indelErrorModel )

// Checks whether anchors are selected outside of repeat regions
BOOST_AUTO_TEST_CASE( test_findingAnchors )
{
    const std::vector<float> knownX = {1, 3};
    const std::vector<float> knownY1 = {1, 3};
    const std::vector<float> knownY2 = {2, 6};
    const std::vector<float> requestX = {2, 4};
    const std::vector<float> requestY1 = {2, 4};
    const std::vector<float> requestY2 = {4, 8};

    for (unsigned ix = 0; ix < requestX.size(); ix++)
    {
        BOOST_CHECK_EQUAL(AdaptiveIndelErrorModel::linearFit(requestX[ix], knownX[0], knownY1[0], knownX[1], knownY1[1]), requestY1[ix]);
        BOOST_CHECK_EQUAL(AdaptiveIndelErrorModel::linearFit(requestX[ix], knownX[0], knownY2[0], knownX[1], knownY2[1]), requestY2[ix]);
    }

}

BOOST_AUTO_TEST_SUITE_END()
