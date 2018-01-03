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

#include "ploidy_util.hh"


BOOST_AUTO_TEST_SUITE( ploidy_util_test )


BOOST_AUTO_TEST_CASE( test_vcf_ploidy_parse )
{
    static const unsigned expectedSampleCount(3);
    static const char* testVcfRecord = "X\t19\t.\t.\t<CNV>\t.\t.\tGRANDMA=COOKIES;END=100;MQ=60\tGT:GQ:CN:FT\t.:.:1:.\t.:.:0:.\t.:.:.:.";

    known_pos_range2 range;
    std::vector<unsigned> ploidy;
    parsePloidyFromVcf(expectedSampleCount, testVcfRecord, range, ploidy);

    BOOST_REQUIRE_EQUAL(range, known_pos_range2(19,100));
    BOOST_REQUIRE_EQUAL(ploidy.size(), 3u);
    BOOST_REQUIRE_EQUAL(ploidy[0], 1u);
    BOOST_REQUIRE_EQUAL(ploidy[1], 0u);
    BOOST_REQUIRE_EQUAL(ploidy[2], 2u);
}


BOOST_AUTO_TEST_SUITE_END()
