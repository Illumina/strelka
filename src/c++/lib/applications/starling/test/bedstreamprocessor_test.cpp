// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
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

#include "test_config.h"


#include "boost/test/unit_test.hpp"
#include "bedstreamprocessor.hh"

BOOST_AUTO_TEST_SUITE( bed_stream_processor_test )

class dummy_variant_sink : public variant_pipe_stage_base
{
public:
    dummy_variant_sink() : variant_pipe_stage_base() {}
    std::vector<std::unique_ptr<digt_site_info>> the_sites;
    std::vector<std::unique_ptr<digt_indel_info>> the_indels;
    void process(std::unique_ptr<site_info> si) override
    {
        the_sites.push_back(downcast<digt_site_info>(std::move(si)));
    }
    void process(std::unique_ptr<indel_info> ii) override
    {
        the_indels.push_back(downcast<digt_indel_info>(std::move(ii)));
    }
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
    indel_key ik;
    ik.pos=50;
    site.reset(new digt_indel_info(ik,
                                   indel_data(ik),
                                   starling_diploid_indel_core(),
                                   starling_indel_report_info(),
                                   starling_indel_sample_report_info()));

    bsp.process(std::move(site));

    BOOST_REQUIRE(!next->the_indels.back()->first().filters.test(VCF_FILTERS::OffTarget));

    ik.pos=105;
    site.reset(new digt_indel_info(ik,
                                   indel_data(ik),
                                   starling_diploid_indel_core(),
                                   starling_indel_report_info(),
                                   starling_indel_sample_report_info()));


    bsp.process(std::move(site));
    BOOST_REQUIRE(next->the_indels.back()->first().filters.test(VCF_FILTERS::OffTarget));

    ik.pos=150;
    site.reset(new digt_indel_info(ik,
                                   indel_data(ik),
                                   starling_diploid_indel_core(),
                                   starling_indel_report_info(),
                                   starling_indel_sample_report_info()));


    bsp.process(std::move(site));
    BOOST_REQUIRE(!next->the_indels.back()->first().filters.test(VCF_FILTERS::OffTarget));

    ik.pos=250;
    site.reset(new digt_indel_info(ik,
                                   indel_data(ik),
                                   starling_diploid_indel_core(),
                                   starling_indel_report_info(),
                                   starling_indel_sample_report_info()));

    bsp.process(std::move(site));
    BOOST_REQUIRE(next->the_indels.back()->first().filters.test(VCF_FILTERS::OffTarget));

}


BOOST_AUTO_TEST_SUITE_END()

