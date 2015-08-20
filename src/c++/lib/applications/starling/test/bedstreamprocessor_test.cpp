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

class dummy_variant_sink : public variant_pipe_stage_base
{
public:
    dummy_variant_sink() : variant_pipe_stage_base() {}
    std::vector<site_info> the_sites;
    std::vector<indel_info> the_indels;
    void process(site_info& si) override
    {
        the_sites.push_back(si);
    }
    void process(indel_info& ii) override
    {
        the_indels.push_back(ii);
    }
};



BOOST_AUTO_TEST_CASE( filters_snps_before_and_after_range )
{
    dummy_variant_sink next;
    bed_stream_processor bsp(TEST_DATA_PATH "/bed_stream_test.bed.gz", "chr1", next);


    site_info site;

    site.pos = 50;
    bsp.process(site);

    BOOST_REQUIRE(!next.the_sites.back().smod.filters.test(VCF_FILTERS::OffTarget));
    site.pos = 105;
    site.smod.clear();

    bsp.process(site);
    BOOST_REQUIRE(next.the_sites.back().smod.filters.test(VCF_FILTERS::OffTarget));

    site.pos = 150;
    site.smod.clear();

    bsp.process(site);
    BOOST_REQUIRE(!next.the_sites.back().smod.filters.test(VCF_FILTERS::OffTarget));
    site.pos = 250;
    site.smod.clear();

    bsp.process(site);
    BOOST_REQUIRE(next.the_sites.back().smod.filters.test(VCF_FILTERS::OffTarget));

}

BOOST_AUTO_TEST_CASE( filters_indels_before_and_after_range )
{
    dummy_variant_sink next;
    bed_stream_processor bsp(TEST_DATA_PATH "/bed_stream_test.bed.gz", "chr1", next);

    indel_info site(20, indel_key(),
                    starling_diploid_indel_core(),
                    starling_indel_report_info(),
                    starling_indel_sample_report_info());

    site.pos = 50;
    bsp.process(site);

    BOOST_REQUIRE(!next.the_indels.back().imod().filters.test(VCF_FILTERS::OffTarget));
    site.pos = 105;
    site.imod().clear();

    bsp.process(site);
    BOOST_REQUIRE(next.the_indels.back().imod().filters.test(VCF_FILTERS::OffTarget));

    site.pos = 150;
    site.imod().clear();

    bsp.process(site);
    BOOST_REQUIRE(!next.the_indels.back().imod().filters.test(VCF_FILTERS::OffTarget));
    site.pos = 250;
    site.imod().clear();

    bsp.process(site);
    BOOST_REQUIRE(next.the_indels.back().imod().filters.test(VCF_FILTERS::OffTarget));

}


BOOST_AUTO_TEST_SUITE_END()

