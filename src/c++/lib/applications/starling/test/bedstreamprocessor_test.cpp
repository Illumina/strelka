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
    std::vector<std::unique_ptr<GermlineDiploidSiteLocusInfo>> the_sites;
    std::vector<std::unique_ptr<GermlineDiploidIndelLocusInfo>> the_indels;
    void process(std::unique_ptr<GermlineSiteLocusInfo> si) override
    {
        the_sites.push_back(downcast<GermlineDiploidSiteLocusInfo>(std::move(si)));
    }
    void process(std::unique_ptr<GermlineIndelLocusInfo> ii) override
    {
        the_indels.push_back(downcast<GermlineDiploidIndelLocusInfo>(std::move(ii)));
    }
};



BOOST_AUTO_TEST_CASE( filters_snps_before_and_after_range )
{
    std::shared_ptr<dummy_variant_sink> next(new dummy_variant_sink);
    bed_stream_processor bsp(TEST_DATA_PATH "/bed_stream_test.bed.gz", "chr1", std::dynamic_pointer_cast<variant_pipe_stage_base>(next));
    const gvcf_options gvcfOptions = gvcf_options();
    const std::string& chromName = "dummy";
    const gvcf_deriv_options gvcfDerivOptions(gvcfOptions, chromName, false);

    auto getNewSite = [&](const pos_t pos)
    {
        std::unique_ptr<GermlineDiploidSiteLocusInfo> site(new GermlineDiploidSiteLocusInfo(gvcfDerivOptions));
        site->pos = pos;
        return site;
    };

    auto testOffTarget = [&]() {
        return next->the_sites.back()->filters.test(GERMLINE_VARIANT_VCF_FILTERS::OffTarget);
    };

    bsp.process(getNewSite(50));
    BOOST_REQUIRE(! testOffTarget());

    bsp.process(getNewSite(105));
    BOOST_REQUIRE(testOffTarget());

    bsp.process(getNewSite(150));
    BOOST_REQUIRE(! testOffTarget());

    bsp.process(getNewSite(250));
    BOOST_REQUIRE(testOffTarget());
}

BOOST_AUTO_TEST_CASE( filters_indels_before_and_after_range )
{
    std::shared_ptr<dummy_variant_sink> next(new dummy_variant_sink);
    bed_stream_processor bsp(TEST_DATA_PATH "/bed_stream_test.bed.gz", "chr1", std::dynamic_pointer_cast<variant_pipe_stage_base>(next));
    const gvcf_options gvcfOptions = gvcf_options();
    const std::string& chromName = "dummy";
    const gvcf_deriv_options gvcfDerivOptions(gvcfOptions, chromName, false);

    auto getNewIndel = [&](const pos_t pos)
    {
        const IndelKey indelKey(pos);
        std::unique_ptr<GermlineDiploidIndelLocusInfo>
        indel(new GermlineDiploidIndelLocusInfo(
                  gvcfDerivOptions,
                  indelKey,
                  IndelData(1,indelKey),
                  GermlineDiploidIndelSimpleGenotypeInfoCore(),
                  starling_indel_report_info(),
                  starling_indel_sample_report_info()));

        return indel;
    };

    auto testOffTarget = [&]() {
        return next->the_indels.back()->filters.test(GERMLINE_VARIANT_VCF_FILTERS::OffTarget);
    };

    bsp.process(getNewIndel(50));
    BOOST_REQUIRE(! testOffTarget());

    bsp.process(getNewIndel(105));
    BOOST_REQUIRE(testOffTarget());

    bsp.process(getNewIndel(150));
    BOOST_REQUIRE(! testOffTarget());

    bsp.process(getNewIndel(250));
    BOOST_REQUIRE(testOffTarget());
}


BOOST_AUTO_TEST_SUITE_END()
