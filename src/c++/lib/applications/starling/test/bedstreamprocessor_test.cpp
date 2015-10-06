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
    std::vector<std::unique_ptr<digt_site_info>> the_sites;
    std::vector<std::unique_ptr<digt_indel_info>> the_indels;
    void process(std::unique_ptr<site_info> si) override { the_sites.push_back(downcast<digt_site_info>(std::move(si))); }
    void process(std::unique_ptr<indel_info> ii) override { the_indels.push_back(downcast<digt_indel_info>(std::move(ii))); }
};



BOOST_AUTO_TEST_CASE( filters_snps_before_and_after_range )
{
    std::shared_ptr<dummy_variant_sink> next(new dummy_variant_sink);
    bed_stream_processor bsp(TEST_DATA_PATH "/bed_stream_test.bed.gz", "chr1", std::dynamic_pointer_cast<variant_pipe_stage_base>(next));


    std::unique_ptr<digt_site_info> site(new digt_site_info());

    site->pos = 50;
    bsp.process(std::move(site));

    BOOST_REQUIRE(!next->the_sites.back()->smod.filters.test(VCF_FILTERS::OffTarget));
    site.reset(new digt_site_info());
    site->pos = 105;

    bsp.process(std::move(site));
    BOOST_REQUIRE(next->the_sites.back()->smod.filters.test(VCF_FILTERS::OffTarget));

    site.reset(new digt_site_info());
    site->pos = 150;

    bsp.process(std::move(site));
    BOOST_REQUIRE(!next->the_sites.back()->smod.filters.test(VCF_FILTERS::OffTarget));

    site.reset(new digt_site_info());
    site->pos = 250;

    bsp.process(std::move(site));
    BOOST_REQUIRE(next->the_sites.back()->smod.filters.test(VCF_FILTERS::OffTarget));

}

BOOST_AUTO_TEST_CASE( filters_indels_before_and_after_range )
{
    std::shared_ptr<dummy_variant_sink> next(new dummy_variant_sink);
    bed_stream_processor bsp(TEST_DATA_PATH "/bed_stream_test.bed.gz", "chr1", std::dynamic_pointer_cast<variant_pipe_stage_base>(next));


    std::unique_ptr<digt_indel_info> site;
    site.reset(new digt_indel_info(50, indel_key(),
         indel_data(indel_key()),
         starling_diploid_indel_core(),
         starling_indel_report_info(),
         starling_indel_sample_report_info()));

    bsp.process(std::move(site));

    BOOST_REQUIRE(!next->the_indels.back()->first().filters.test(VCF_FILTERS::OffTarget));

    site.reset(new digt_indel_info(105, indel_key(),
             indel_data(indel_key()),
             starling_diploid_indel_core(),
             starling_indel_report_info(),
             starling_indel_sample_report_info()));


    bsp.process(std::move(site));
    BOOST_REQUIRE(next->the_indels.back()->first().filters.test(VCF_FILTERS::OffTarget));

    site.reset(new digt_indel_info(150, indel_key(),
            indel_data(indel_key()),
             starling_diploid_indel_core(),
             starling_indel_report_info(),
             starling_indel_sample_report_info()));


    bsp.process(std::move(site));
    BOOST_REQUIRE(!next->the_indels.back()->first().filters.test(VCF_FILTERS::OffTarget));
    site.reset(new digt_indel_info(250, indel_key(),
            indel_data(indel_key()),
             starling_diploid_indel_core(),
             starling_indel_report_info(),
             starling_indel_sample_report_info()));

    bsp.process(std::move(site));
    BOOST_REQUIRE(next->the_indels.back()->first().filters.test(VCF_FILTERS::OffTarget));

}


BOOST_AUTO_TEST_SUITE_END()

