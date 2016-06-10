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
    std::vector<std::unique_ptr<GermlineDiploidSiteCallInfo>> the_sites;
    std::vector<std::unique_ptr<GermlineDiploidIndelCallInfo>> the_indels;
    void process(std::unique_ptr<GermlineSiteCallInfo> si) override
    {
        the_sites.push_back(downcast<GermlineDiploidSiteCallInfo>(std::move(si)));
    }
    void process(std::unique_ptr<GermlineIndelCallInfo> ii) override
    {
        the_indels.push_back(downcast<GermlineDiploidIndelCallInfo>(std::move(ii)));
    }
};



BOOST_AUTO_TEST_CASE( filters_snps_before_and_after_range )
{
    std::shared_ptr<dummy_variant_sink> next(new dummy_variant_sink);
    bed_stream_processor bsp(TEST_DATA_PATH "/bed_stream_test.bed.gz", "chr1", std::dynamic_pointer_cast<variant_pipe_stage_base>(next));


    std::unique_ptr<GermlineDiploidSiteCallInfo> site(new GermlineDiploidSiteCallInfo());

    site->pos = 50;
    bsp.process(std::move(site));

    BOOST_REQUIRE(!next->the_sites.back()->smod.filters.test(GERMLINE_VARIANT_VCF_FILTERS::OffTarget));
    site.reset(new GermlineDiploidSiteCallInfo());
    site->pos = 105;

    bsp.process(std::move(site));
    BOOST_REQUIRE(next->the_sites.back()->smod.filters.test(GERMLINE_VARIANT_VCF_FILTERS::OffTarget));

    site.reset(new GermlineDiploidSiteCallInfo());
    site->pos = 150;

    bsp.process(std::move(site));
    BOOST_REQUIRE(!next->the_sites.back()->smod.filters.test(GERMLINE_VARIANT_VCF_FILTERS::OffTarget));

    site.reset(new GermlineDiploidSiteCallInfo());
    site->pos = 250;

    bsp.process(std::move(site));
    BOOST_REQUIRE(next->the_sites.back()->smod.filters.test(GERMLINE_VARIANT_VCF_FILTERS::OffTarget));

}

BOOST_AUTO_TEST_CASE( filters_indels_before_and_after_range )
{
    std::shared_ptr<dummy_variant_sink> next(new dummy_variant_sink);
    bed_stream_processor bsp(TEST_DATA_PATH "/bed_stream_test.bed.gz", "chr1", std::dynamic_pointer_cast<variant_pipe_stage_base>(next));


    std::unique_ptr<GermlineDiploidIndelCallInfo> site;
    IndelKey indelKey;
    indelKey.pos=50;
    site.reset(new GermlineDiploidIndelCallInfo(indelKey,
                                                IndelData(1,indelKey),
                                                GermlineDiploidIndelSimpleGenotypeInfoCore(),
                                                starling_indel_report_info(),
                                                starling_indel_sample_report_info()));

    bsp.process(std::move(site));

    BOOST_REQUIRE(!next->the_indels.back()->first().filters.test(GERMLINE_VARIANT_VCF_FILTERS::OffTarget));

    indelKey.pos=105;
    site.reset(new GermlineDiploidIndelCallInfo(indelKey,
                                                IndelData(1,indelKey),
                                                GermlineDiploidIndelSimpleGenotypeInfoCore(),
                                                starling_indel_report_info(),
                                                starling_indel_sample_report_info()));


    bsp.process(std::move(site));
    BOOST_REQUIRE(next->the_indels.back()->first().filters.test(GERMLINE_VARIANT_VCF_FILTERS::OffTarget));

    indelKey.pos=150;
    site.reset(new GermlineDiploidIndelCallInfo(indelKey,
                                                IndelData(1,indelKey),
                                                GermlineDiploidIndelSimpleGenotypeInfoCore(),
                                                starling_indel_report_info(),
                                                starling_indel_sample_report_info()));


    bsp.process(std::move(site));
    BOOST_REQUIRE(!next->the_indels.back()->first().filters.test(GERMLINE_VARIANT_VCF_FILTERS::OffTarget));

    indelKey.pos=250;
    site.reset(new GermlineDiploidIndelCallInfo(indelKey,
                                                IndelData(1,indelKey),
                                                GermlineDiploidIndelSimpleGenotypeInfoCore(),
                                                starling_indel_report_info(),
                                                starling_indel_sample_report_info()));

    bsp.process(std::move(site));
    BOOST_REQUIRE(next->the_indels.back()->first().filters.test(GERMLINE_VARIANT_VCF_FILTERS::OffTarget));
}


BOOST_AUTO_TEST_SUITE_END()

