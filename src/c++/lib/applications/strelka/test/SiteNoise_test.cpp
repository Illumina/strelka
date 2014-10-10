// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
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

