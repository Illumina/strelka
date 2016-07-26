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
#include "boost/test/unit_test.hpp"
#include "starling_continuous_variant_caller.hh"
#include <vector>

BOOST_AUTO_TEST_SUITE( continuous_variant_caller )

BOOST_AUTO_TEST_CASE( qscore_calculation )
{
    std::vector<std::vector<unsigned>> SampleValues
    {
        //coverage,var calls, expected Q score
        {100,0,0},
        {100,1,2},
        {100,5,24},
        {200,10,43},
        {500,25,98},
        {5000,250,100}
    };

    for (auto& item :  SampleValues)
    {
        auto FinalQValue = starling_continuous_variant_caller::poisson_qscore(item[1], item[0], 20, 100);
        BOOST_REQUIRE_EQUAL(FinalQValue, item[2]);
    }
}

BOOST_AUTO_TEST_CASE( call_from_counts )
{
    starling_base_options opt;
    snp_pos_info pileup;

    for (int i=0; i<10; i++)
        pileup.calls.emplace_back(base_to_id('C'), 30, false, 0, 0, false, false, false, false, false);
    for (int i=0; i<20; i++)
        pileup.calls.emplace_back(base_to_id('G'), 30, false, 0, 0, false, false, false, false, false);
//    for (int i=0;i<30; i++)
//        pileup.calls.emplace_back(base_to_id('T'), 30, false, 0, 0, false, false, false, false, false);
    for (int i=0; i<40; i++)
        pileup.calls.emplace_back(base_to_id('A'), 30, false, 0, 0, false, false, false, false, false);

    const unsigned sampleCount(1);
    GermlineContinuousSiteLocusInfo info(sampleCount, 10, 'A', pileup, opt.min_qscore, opt.min_het_vf);

    starling_continuous_variant_caller::position_snp_call_continuous(opt, pileup, BASE_ID::C, false, info);
    auto C = std::find_if(info.altAlleles.begin(), info.altAlleles.end(), [&](const GermlineContinuousSiteAlleleInfo& call)
    {
        return call._base == BASE_ID::C;
    });
    BOOST_REQUIRE(C != info.altAlleles.end());
    BOOST_REQUIRE_EQUAL(10, C->_alleleDepth);
    BOOST_REQUIRE_EQUAL(40, C->gqx);
    BOOST_REQUIRE_EQUAL(70, C->_totalDepth);

    info.altAlleles.clear();
    starling_continuous_variant_caller::position_snp_call_continuous(opt, pileup, BASE_ID::G, false, info);
    auto G = std::find_if(info.altAlleles.begin(), info.altAlleles.end(), [&](const GermlineContinuousSiteAlleleInfo& call)
    {
        return call._base == BASE_ID::G;
    });
    BOOST_REQUIRE(G != info.altAlleles.end());
    BOOST_REQUIRE_EQUAL(20, G->_alleleDepth);
    BOOST_REQUIRE_EQUAL(40, G->gqx);
    BOOST_REQUIRE_EQUAL(70, G->_totalDepth);

    info.altAlleles.clear();
    starling_continuous_variant_caller::position_snp_call_continuous(opt, pileup, BASE_ID::T, false, info);
    auto T = std::find_if(info.altAlleles.begin(), info.altAlleles.end(), [&](const GermlineContinuousSiteAlleleInfo& call)
    {
        return call._base == BASE_ID::T;
    });
    BOOST_REQUIRE(T == info.altAlleles.end());

    info.altAlleles.clear();
    starling_continuous_variant_caller::position_snp_call_continuous(opt, pileup, BASE_ID::A, true, info);
    auto A = std::find_if(info.altAlleles.begin(), info.altAlleles.end(), [&](const GermlineContinuousSiteAlleleInfo& call)
    {
        return call._base == BASE_ID::A;
    });
    BOOST_REQUIRE(A != info.altAlleles.end());
    BOOST_REQUIRE_EQUAL(40, A->_alleleDepth);
    BOOST_REQUIRE_EQUAL(40, A->gqx);
    BOOST_REQUIRE_EQUAL(70, A->_totalDepth);
}

BOOST_AUTO_TEST_CASE( do_not_call_low_vf )
{
    starling_base_options opt;
    snp_pos_info pileup;

    opt.min_het_vf = 0.03; // 3% threshold

    // insert C at 2%, should not be called
    for (int i=0; i<2; i++)
        pileup.calls.emplace_back(base_to_id('C'), 30, false, 0, 0, false, false, false, false, false);

    for (int i=0; i<5; i++)
        pileup.calls.emplace_back(base_to_id('T'), 30, false, 0, 0, false, false, false, false, false);

    for (int i=0; i<93; i++)
        pileup.calls.emplace_back(base_to_id('A'), 30, false, 0, 0, false, false, false, false, false);

    const unsigned sampleCount(1);
    GermlineContinuousSiteLocusInfo info(sampleCount, 10, 'A', pileup, opt.min_qscore, opt.min_het_vf);

    starling_continuous_variant_caller::position_snp_call_continuous(opt, pileup, BASE_ID::C, false, info);
    auto C = std::find_if(info.altAlleles.begin(), info.altAlleles.end(), [&](const GermlineContinuousSiteAlleleInfo& call)
    {
        return call._base == BASE_ID::C;
    });
    BOOST_REQUIRE(C == info.altAlleles.end());

    info.altAlleles.clear();
    starling_continuous_variant_caller::position_snp_call_continuous(opt, pileup, BASE_ID::A, true, info);
    auto A = std::find_if(info.altAlleles.begin(), info.altAlleles.end(), [&](const GermlineContinuousSiteAlleleInfo& call)
    {
        return call._base == BASE_ID::A;
    });
    BOOST_REQUIRE(A != info.altAlleles.end());
    BOOST_REQUIRE_EQUAL(93, A->_alleleDepth);
    BOOST_REQUIRE_EQUAL(40, A->gqx);
    BOOST_REQUIRE_EQUAL(100, A->_totalDepth);

    info.altAlleles.clear();
    starling_continuous_variant_caller::position_snp_call_continuous(opt, pileup, BASE_ID::T, false, info);
    auto T = std::find_if(info.altAlleles.begin(), info.altAlleles.end(), [&](const GermlineContinuousSiteAlleleInfo& call)
    {
        return call._base == BASE_ID::T;
    });
    BOOST_REQUIRE(T != info.altAlleles.end());
    BOOST_REQUIRE_EQUAL(5, T->_alleleDepth);
    BOOST_REQUIRE_EQUAL(13, T->gqx);
    BOOST_REQUIRE_EQUAL(100, T->_totalDepth);
}

BOOST_AUTO_TEST_CASE( gt_forced_output_calculated_correctly )
{
    starling_base_options opt;
    snp_pos_info pileup;

    opt.min_het_vf = 0.03; // 3% threshold

    // insert C at 2%, should not be called normally, but it will be if forced
    for (int i = 0; i < 2; i++)
        pileup.calls.emplace_back(base_to_id('C'), 30, false, 0, 0, false, false, false, false, false);

    for (int i = 0; i < 5; i++)
        pileup.calls.emplace_back(base_to_id('T'), 30, false, 0, 0, false, false, false, false, false);

    for (int i = 0; i < 93; i++)
        pileup.calls.emplace_back(base_to_id('A'), 30, false, 0, 0, false, false, false, false, false);

    const unsigned sampleCount(1);
    GermlineContinuousSiteLocusInfo info(sampleCount, 10, 'A', pileup, opt.min_qscore, opt.min_het_vf, true);

    starling_continuous_variant_caller::position_snp_call_continuous(opt, pileup, BASE_ID::C, true, info);
    auto C = std::find_if(info.altAlleles.begin(), info.altAlleles.end(),
                          [&](const GermlineContinuousSiteAlleleInfo& call) {
                              return call._base == BASE_ID::C;
                          });
    BOOST_REQUIRE(C != info.altAlleles.end());
    // below min het vf GT should be 0/0 (STAR-66)
    BOOST_REQUIRE(0 == strcmp(info.get_gt(*C), "0/0"));

    info.altAlleles.clear();
    starling_continuous_variant_caller::position_snp_call_continuous(opt, pileup, BASE_ID::A, true, info);
    auto A = std::find_if(info.altAlleles.begin(), info.altAlleles.end(),
                          [&](const GermlineContinuousSiteAlleleInfo& call) {
                              return call._base == BASE_ID::A;
                          });
    BOOST_REQUIRE(A != info.altAlleles.end());
    BOOST_REQUIRE_EQUAL(93, A->_alleleDepth);
    BOOST_REQUIRE_EQUAL(40, A->gqx);
    BOOST_REQUIRE_EQUAL(100, A->_totalDepth);
    BOOST_REQUIRE(0 == strcmp(info.get_gt(*A), "0/0"));

    info.altAlleles.clear();
    starling_continuous_variant_caller::position_snp_call_continuous(opt, pileup, BASE_ID::T, true, info);
    auto T = std::find_if(info.altAlleles.begin(), info.altAlleles.end(),
                          [&](const GermlineContinuousSiteAlleleInfo& call) {
                              return call._base == BASE_ID::T;
                          });
    BOOST_REQUIRE(T != info.altAlleles.end());
    BOOST_REQUIRE_EQUAL(5, T->_alleleDepth);
    BOOST_REQUIRE_EQUAL(13, T->gqx);
    BOOST_REQUIRE_EQUAL(100, T->_totalDepth);
    BOOST_REQUIRE(0 == strcmp(info.get_gt(*T), "0/1"));

    // G will be called since it is forced output
    // below min het vf GT should be 0/0 (STAR-66)
    info.altAlleles.clear();
    starling_continuous_variant_caller::position_snp_call_continuous(opt, pileup, BASE_ID::G, true, info);
    auto G = std::find_if(info.altAlleles.begin(), info.altAlleles.end(),
                          [&](const GermlineContinuousSiteAlleleInfo& call) {
                              return call._base == BASE_ID::G;
                          });
    BOOST_REQUIRE(G != info.altAlleles.end());
    BOOST_REQUIRE_EQUAL(0, G->_alleleDepth);
    BOOST_REQUIRE_EQUAL(100, G->_totalDepth);
    BOOST_REQUIRE(0 == strcmp(info.get_gt(*G), "0/0"));
}


BOOST_AUTO_TEST_CASE( homalt_called_correctly )
{
    starling_base_options opt;
    snp_pos_info pileup;

    opt.min_het_vf = 0.03; // 3% threshold

    // insert A at 2%, should not be called normally
    for (int i=0; i<2; i++)
        pileup.calls.emplace_back(base_to_id('A'), 30, false, 0, 0, false, false, false, false, false);

    for (int i=0; i<98; i++)
        pileup.calls.emplace_back(base_to_id('T'), 30, false, 0, 0, false, false, false, false, false);

    const unsigned sampleCount(1);
    GermlineContinuousSiteLocusInfo info(sampleCount, 10, 'A', pileup, opt.min_qscore, opt.min_het_vf);

    starling_continuous_variant_caller::position_snp_call_continuous(opt, pileup, BASE_ID::T, false, info);
    auto T = std::find_if(info.altAlleles.begin(), info.altAlleles.end(), [&](const GermlineContinuousSiteAlleleInfo& call)
    {
        return call._base == BASE_ID::T;
    });
    BOOST_REQUIRE(T != info.altAlleles.end());
    BOOST_REQUIRE_EQUAL(98, T->_alleleDepth);
    BOOST_REQUIRE_EQUAL(100, T->_totalDepth);
    BOOST_REQUIRE(0 == strcmp(info.get_gt(*T), "1/1"));

    BOOST_REQUIRE_EQUAL(1, info.altAlleles.size());
}


BOOST_AUTO_TEST_CASE( compute_sb )
{
    // TODO: diploid one returns -11776.1. chrM 14213
    BOOST_REQUIRE_CLOSE_FRACTION(0.246095, starling_continuous_variant_caller::strand_bias(4769,4058, 2, 4, .01), 0.0001);
    // chr1 16892310, SB: -65.243244, fwdAlt: 0, revAlt:62, fwdOther:190, revOther:8
    // diplod one returns 131.7
    BOOST_REQUIRE_CLOSE_FRACTION(-65.243244, starling_continuous_variant_caller::strand_bias(0,62, 190, 8, .01), 0.0001);
}


BOOST_AUTO_TEST_SUITE_END()
