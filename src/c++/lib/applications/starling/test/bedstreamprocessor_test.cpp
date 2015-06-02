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

#include "test_config.h"


#include "boost/test/unit_test.hpp"
#include "bedstreamprocessor.cpp"

BOOST_AUTO_TEST_SUITE( bed_stream_processor_test )



BOOST_AUTO_TEST_CASE( filters_snps_before_and_after_range )
{
    bed_stream_processor bsp(TEST_DATA_PATH "/bed_stream_test.bed.gz", "chr1");

    site_info site;

    site.pos = 50;
    bsp.process(site);

    BOOST_REQUIRE(!site.smod.filters.test(VCF_FILTERS::OffTarget));
    site.pos = 105;
    site.smod.clear();

    bsp.process(site);
    BOOST_REQUIRE(site.smod.filters.test(VCF_FILTERS::OffTarget));

    site.pos = 150;
    site.smod.clear();

    bsp.process(site);
    BOOST_REQUIRE(!site.smod.filters.test(VCF_FILTERS::OffTarget));
    site.pos = 250;
    site.smod.clear();

    bsp.process(site);
    BOOST_REQUIRE(site.smod.filters.test(VCF_FILTERS::OffTarget));

}

BOOST_AUTO_TEST_CASE( filters_indels_before_and_after_range )
{
    bed_stream_processor bsp(TEST_DATA_PATH "/bed_stream_test.bed.gz", "chr1");

    indel_info site;

    site.pos = 50;
    bsp.process(site);

    BOOST_REQUIRE(!site.imod.filters.test(VCF_FILTERS::OffTarget));
    site.pos = 105;
    site.imod.clear();

    bsp.process(site);
    BOOST_REQUIRE(site.imod.filters.test(VCF_FILTERS::OffTarget));

    site.pos = 150;
    site.imod.clear();

    bsp.process(site);
    BOOST_REQUIRE(!site.imod.filters.test(VCF_FILTERS::OffTarget));
    site.pos = 250;
    site.imod.clear();

    bsp.process(site);
    BOOST_REQUIRE(site.imod.filters.test(VCF_FILTERS::OffTarget));

}


BOOST_AUTO_TEST_SUITE_END()

