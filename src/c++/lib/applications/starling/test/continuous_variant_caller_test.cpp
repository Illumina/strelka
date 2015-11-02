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
#include "starling_continuous_variant_caller.cpp"
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
        BOOST_CHECK_EQUAL(FinalQValue, item[2]);
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


    continuous_site_info info(10, 'A', pileup, opt.min_qscore, opt.min_het_vf);

    starling_continuous_variant_caller::position_snp_call_continuous(opt, pileup, info);
    auto C = std::find_if(info.calls.begin(), info.calls.end(), [&](const continuous_site_call& call)
    {
        return call._base == BASE_ID::C;
    });
    BOOST_CHECK(C != info.calls.end());
    BOOST_CHECK_EQUAL(10, C->_alleleDepth);
    BOOST_CHECK_EQUAL(40, C->gq);
    BOOST_CHECK_EQUAL(70, C->_totalDepth);

    auto G = std::find_if(info.calls.begin(), info.calls.end(), [&](const continuous_site_call& call)
    {
        return call._base == BASE_ID::G;
    });
    BOOST_CHECK(G != info.calls.end());
    BOOST_CHECK_EQUAL(20, G->_alleleDepth);
    BOOST_CHECK_EQUAL(40, G->gq);
    BOOST_CHECK_EQUAL(70, G->_totalDepth);

    auto T = std::find_if(info.calls.begin(), info.calls.end(), [&](const continuous_site_call& call)
    {
        return call._base == BASE_ID::T;
    });
    BOOST_CHECK(T == info.calls.end());

    auto A = std::find_if(info.calls.begin(), info.calls.end(), [&](const continuous_site_call& call)
    {
        return call._base == BASE_ID::A;
    });
    BOOST_CHECK(A != info.calls.end());
    BOOST_CHECK_EQUAL(40, A->_alleleDepth);
    BOOST_CHECK_EQUAL(40, A->gq);
    BOOST_CHECK_EQUAL(70, A->_totalDepth);
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


    continuous_site_info info(10, 'A', pileup, opt.min_qscore, opt.min_het_vf);

    starling_continuous_variant_caller::position_snp_call_continuous(opt, pileup, info);
    auto C = std::find_if(info.calls.begin(), info.calls.end(), [&](const continuous_site_call& call)
    {
        return call._base == BASE_ID::C;
    });
    BOOST_CHECK(C == info.calls.end());

    auto A = std::find_if(info.calls.begin(), info.calls.end(), [&](const continuous_site_call& call)
    {
        return call._base == BASE_ID::A;
    });
    BOOST_CHECK(A != info.calls.end());
    BOOST_CHECK_EQUAL(93, A->_alleleDepth);
    BOOST_CHECK_EQUAL(40, A->gq);
    BOOST_CHECK_EQUAL(100, A->_totalDepth);

    auto T = std::find_if(info.calls.begin(), info.calls.end(), [&](const continuous_site_call& call)
    {
        return call._base == BASE_ID::T;
    });
    BOOST_CHECK(T != info.calls.end());
    BOOST_CHECK_EQUAL(5, T->_alleleDepth);
    BOOST_CHECK_EQUAL(13, T->gq);
    BOOST_CHECK_EQUAL(100, T->_totalDepth);

    BOOST_CHECK_EQUAL(2, info.calls.size());
}

BOOST_AUTO_TEST_CASE( gt_forced_output_calculated_correctly )
{
    starling_base_options opt;
    snp_pos_info pileup;

    opt.min_het_vf = 0.03; // 3% threshold

    // insert C at 2%, should not be called normally, but it will be if forced
    for (int i=0; i<2; i++)
        pileup.calls.emplace_back(base_to_id('C'), 30, false, 0, 0, false, false, false, false, false);

    for (int i=0; i<5; i++)
        pileup.calls.emplace_back(base_to_id('T'), 30, false, 0, 0, false, false, false, false, false);

    for (int i=0; i<93; i++)
        pileup.calls.emplace_back(base_to_id('A'), 30, false, 0, 0, false, false, false, false, false);


    continuous_site_info info(10, 'A', pileup, opt.min_qscore, opt.min_het_vf, true);

    starling_continuous_variant_caller::position_snp_call_continuous(opt, pileup, info);
    auto C = std::find_if(info.calls.begin(), info.calls.end(), [&](const continuous_site_call& call)
    {
        return call._base == BASE_ID::C;
    });
    BOOST_CHECK(C != info.calls.end());
    BOOST_CHECK(0 == strcmp(info.get_gt(*C), "0/1"));

    auto A = std::find_if(info.calls.begin(), info.calls.end(), [&](const continuous_site_call& call)
    {
        return call._base == BASE_ID::A;
    });
    BOOST_CHECK(A != info.calls.end());
    BOOST_CHECK_EQUAL(93, A->_alleleDepth);
    BOOST_CHECK_EQUAL(40, A->gq);
    BOOST_CHECK_EQUAL(100, A->_totalDepth);
    BOOST_CHECK(0 == strcmp(info.get_gt(*A), "0/0"));

    auto T = std::find_if(info.calls.begin(), info.calls.end(), [&](const continuous_site_call& call)
    {
        return call._base == BASE_ID::T;
    });
    BOOST_CHECK(T != info.calls.end());
    BOOST_CHECK_EQUAL(5, T->_alleleDepth);
    BOOST_CHECK_EQUAL(13, T->gq);
    BOOST_CHECK_EQUAL(100, T->_totalDepth);
    BOOST_CHECK(0 == strcmp(info.get_gt(*T), "0/1"));

    // G will be called since it is forced output
    auto G = std::find_if(info.calls.begin(), info.calls.end(), [&](const continuous_site_call& call)
    {
        return call._base == BASE_ID::G;
    });
    BOOST_CHECK(G != info.calls.end());
    BOOST_CHECK_EQUAL(0, G->_alleleDepth);
    BOOST_CHECK_EQUAL(100, G->_totalDepth);
    BOOST_CHECK(0 == strcmp(info.get_gt(*G), "0/1"));

    BOOST_CHECK_EQUAL(4, info.calls.size());
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


    continuous_site_info info(10, 'A', pileup, opt.min_qscore, opt.min_het_vf);

    starling_continuous_variant_caller::position_snp_call_continuous(opt, pileup, info);

    auto T = std::find_if(info.calls.begin(), info.calls.end(), [&](const continuous_site_call& call)
    {
        return call._base == BASE_ID::T;
    });
    BOOST_CHECK(T != info.calls.end());
    BOOST_CHECK_EQUAL(98, T->_alleleDepth);
    BOOST_CHECK_EQUAL(100, T->_totalDepth);
    BOOST_CHECK(0 == strcmp(info.get_gt(*T), "1/1"));

    BOOST_CHECK_EQUAL(1, info.calls.size());
}


BOOST_AUTO_TEST_CASE( compute_sb )
{
    // TODO: diploid one returns -11776.1. chrM 14213
    BOOST_CHECK_CLOSE_FRACTION(0.246095, starling_continuous_variant_caller::strand_bias(4769,4058, 2, 4, .01), 0.0001);
    // chr1 16892310, SB: -65.243244, fwdAlt: 0, revAlt:62, fwdOther:190, revOther:8
    // diplod one returns 131.7
    BOOST_CHECK_CLOSE_FRACTION(-65.243244, starling_continuous_variant_caller::strand_bias(0,62, 190, 8, .01), 0.0001);



}



BOOST_AUTO_TEST_SUITE_END()
