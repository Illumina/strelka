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

#include "SiteNoise.hh"


BOOST_AUTO_TEST_SUITE( SiteNoise_test )


BOOST_AUTO_TEST_CASE( test_SiteNoise )
{
    static const char line[] = "chr1\t1941613\t.\tT\tA\t.\t.\t.\tDP:AD\t.:.\t31:30,1\t34:33,3\t.:.\t.:.\t.:.\t.:.\t.:.\t.:.\n";

    SiteNoise sn;
    set_noise_from_vcf(line,sn);

    BOOST_REQUIRE_EQUAL(sn.total,9);
    BOOST_REQUIRE_EQUAL(sn.noise,2);
    BOOST_REQUIRE_EQUAL(sn.noise2,1);
}


BOOST_AUTO_TEST_SUITE_END()

