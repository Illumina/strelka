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

#include "starling_read_align.cpp"



BOOST_AUTO_TEST_SUITE( starling_read_align )

// test results of placing a single indel into a standard background:
static
candidate_alignment
test_indel_placement(
    const IndelKey& indelKey,
    const pos_t read_start_pos = 0)
{
    const pos_t ref_start_pos(1000);
    const bool is_fwd_strand(true);
    const unsigned read_length(100);

    indel_set_t iset;
    const IndelKey ikfixed(1075,INDEL::INDEL,1);
    iset.insert(ikfixed);
    iset.insert(indelKey);

    return make_start_pos_alignment(ref_start_pos,read_start_pos,is_fwd_strand,read_length,iset);
}


static
void
path_compare(const std::string& expect,
             const ALIGNPATH::path_t& result)
{

    using namespace ALIGNPATH;

    path_t tmp_path;
    cigar_to_apath(expect.c_str(),tmp_path);
    BOOST_CHECK_EQUAL(tmp_path,result);
}


BOOST_AUTO_TEST_CASE( test_make_start_pos_alignment )
{
    const char* insertSeq10 = "AAAAACCCCC";
    {
        // basic delete
        IndelKey ik(1050,INDEL::INDEL,10);
        candidate_alignment cal = test_indel_placement(ik);
        path_compare("50M10D15M1D35M",cal.al.path);
    }

    {
        // basic insert
        IndelKey ik(1050,INDEL::INDEL, 0, insertSeq10);
        candidate_alignment cal = test_indel_placement(ik);
        path_compare("50M10I25M1D15M",cal.al.path);
    }

    {
        // basic swap
        IndelKey ik(1050,INDEL::INDEL, 5, insertSeq10);
        candidate_alignment cal = test_indel_placement(ik);
        path_compare("50M5D10I20M1D20M",cal.al.path);
    }


    {
        // trailing edge insert
        IndelKey ik(1091,INDEL::INDEL, 0, insertSeq10);
        candidate_alignment cal = test_indel_placement(ik);
        path_compare("75M1D15M10I",cal.al.path);
        BOOST_CHECK_EQUAL(cal.trailing_indel_key,ik);
    }

    {
        // trailing edge insert
        IndelKey ik(1096,INDEL::INDEL, 0, insertSeq10);
        candidate_alignment cal = test_indel_placement(ik);
        path_compare("75M1D20M5I",cal.al.path);
        BOOST_CHECK_EQUAL(cal.trailing_indel_key,ik);
    }

    {
        // trailing edge insert miss
        IndelKey ik(1101,INDEL::INDEL, 0, insertSeq10);
        candidate_alignment cal = test_indel_placement(ik);
        path_compare("75M1D25M",cal.al.path);
        BOOST_CHECK_EQUAL(cal.trailing_indel_key.type,INDEL::NONE);
    }

    {
        // leading edge insert miss
        IndelKey ik(1000,INDEL::INDEL, 0, insertSeq10);
        candidate_alignment cal = test_indel_placement(ik);
        path_compare("75M1D25M",cal.al.path);
        BOOST_CHECK_EQUAL(cal.leading_indel_key.type,INDEL::NONE);
        BOOST_CHECK_EQUAL(cal.al.pos,1000);
    }

    {
        // leading edge insert
        IndelKey ik(1000,INDEL::INDEL, 0, insertSeq10);
        candidate_alignment cal = test_indel_placement(ik,5);
        path_compare("5I75M1D20M",cal.al.path);
        BOOST_CHECK_EQUAL(cal.leading_indel_key,ik);
        BOOST_CHECK_EQUAL(cal.al.pos,1000);
    }

    {
        // leading edge swap
        IndelKey ik(1000,INDEL::INDEL, 5, insertSeq10);
        candidate_alignment cal = test_indel_placement(ik,5);
        path_compare("5I5D70M1D25M",cal.al.path);
        BOOST_CHECK_EQUAL(cal.leading_indel_key,ik);
        BOOST_CHECK_EQUAL(cal.al.pos,1000);
    }

    {
        // leading edge delete:
        IndelKey ik(1000,INDEL::INDEL,10);
        candidate_alignment cal = test_indel_placement(ik);
        path_compare("10D65M1D35M",cal.al.path);
        BOOST_CHECK_EQUAL(cal.leading_indel_key,ik);
        BOOST_CHECK_EQUAL(cal.al.pos,1000);
    }

    {
        // trailing edge delete:
        IndelKey ik(1101,INDEL::INDEL,10);
        candidate_alignment cal = test_indel_placement(ik);
        path_compare("75M1D25M10D",cal.al.path);
        BOOST_CHECK_EQUAL(cal.trailing_indel_key,ik);
        BOOST_CHECK_EQUAL(cal.al.pos,1000);
    }

    {
        // trailing off-edge delete:
        IndelKey ik(1102,INDEL::INDEL,10);
        candidate_alignment cal = test_indel_placement(ik);
        path_compare("75M1D25M",cal.al.path);
        BOOST_CHECK_EQUAL(cal.leading_indel_key.type,INDEL::NONE);
        BOOST_CHECK_EQUAL(cal.al.pos,1000);
    }
}



/// test results of placing a single indel into a standard background:
///
/// return (ref_start_pos,read_start_pos)
///
static
std::pair<int,int>
test_end_pin_indel_placement(
    const IndelKey& indelKey,
    const pos_t read_end_pos = 100)
{
    const pos_t ref_end_pos(1100);
    const unsigned read_length(100);

    indel_set_t iset;
    const IndelKey ikfixed(1075,INDEL::INDEL,1);
    iset.insert(ikfixed);
    iset.insert(indelKey);

    pos_t ref_start_pos(0);
    pos_t read_start_pos(0);
    get_end_pin_start_pos(iset,read_length,ref_end_pos,read_end_pos,ref_start_pos,read_start_pos);

    return std::make_pair(static_cast<int>(ref_start_pos),
                          static_cast<int>(read_start_pos));
}



BOOST_AUTO_TEST_CASE( test_end_pin_start_pos )
{
    const char* insertSeq10 = "AAAAACCCCC";
    {
        // basic delete:
        IndelKey ik(1050,INDEL::INDEL,10);

        const std::pair<int,int> res = test_end_pin_indel_placement(ik);
        BOOST_CHECK_EQUAL(res.first,989);
        BOOST_CHECK_EQUAL(res.second,0);
    }

    {
        // basic insert:
        IndelKey ik(1050,INDEL::INDEL, 0, insertSeq10);

        const std::pair<int,int> res = test_end_pin_indel_placement(ik);
        BOOST_CHECK_EQUAL(res.first,1009);
        BOOST_CHECK_EQUAL(res.second,0);
    }

    {
        // basic swap:
        IndelKey ik(1050,INDEL::INDEL, 5, insertSeq10);

        const std::pair<int,int> res = test_end_pin_indel_placement(ik);
        BOOST_CHECK_EQUAL(res.first,1004);
        BOOST_CHECK_EQUAL(res.second,0);
    }

    {
        // leading edge insert:
        IndelKey ik(1005,INDEL::INDEL, 0, insertSeq10);

        const std::pair<int,int> res = test_end_pin_indel_placement(ik);
        BOOST_CHECK_EQUAL(res.first,1005);
        BOOST_CHECK_EQUAL(res.second,6);
    }

    {
        // leading edge insert miss:
        IndelKey ik(999,INDEL::INDEL, 0, insertSeq10);

        const std::pair<int,int> res = test_end_pin_indel_placement(ik);
        BOOST_CHECK_EQUAL(res.first,999);
        BOOST_CHECK_EQUAL(res.second,0);
    }

    {
        // leading edge insert far miss:
        IndelKey ik(99,INDEL::INDEL, 0, insertSeq10);

        const std::pair<int,int> res = test_end_pin_indel_placement(ik);
        BOOST_CHECK_EQUAL(res.first,999);
        BOOST_CHECK_EQUAL(res.second,0);
    }

    {
        // leading edge swap:
        IndelKey ik(1000,INDEL::INDEL, 5, insertSeq10);

        const std::pair<int,int> res = test_end_pin_indel_placement(ik);
        BOOST_CHECK_EQUAL(res.first,1000);
        BOOST_CHECK_EQUAL(res.second,6);
    }

    {
        // trailing edge insert:
        IndelKey ik(1100,INDEL::INDEL, 0, insertSeq10);

        const std::pair<int,int> res = test_end_pin_indel_placement(ik,95);
        BOOST_CHECK_EQUAL(res.first,1004);
        BOOST_CHECK_EQUAL(res.second,0);
    }

    {
        // trailing edge insert edge miss:
        IndelKey ik(1100,INDEL::INDEL, 0, insertSeq10);

        const std::pair<int,int> res = test_end_pin_indel_placement(ik);
        BOOST_CHECK_EQUAL(res.first,999);
        BOOST_CHECK_EQUAL(res.second,0);
    }

    {
        // trail edge insert miss:
        IndelKey ik(1110,INDEL::INDEL, 0, insertSeq10);

        const std::pair<int,int> res = test_end_pin_indel_placement(ik);
        BOOST_CHECK_EQUAL(res.first,999);
        BOOST_CHECK_EQUAL(res.second,0);
    }

    {
        // trailing edge insert far miss:
        IndelKey ik(2000,INDEL::INDEL, 0, insertSeq10);

        const std::pair<int,int> res = test_end_pin_indel_placement(ik);
        BOOST_CHECK_EQUAL(res.first,999);
        BOOST_CHECK_EQUAL(res.second,0);
    }

    {
        // near (but not over) trailing edge swap:
        IndelKey ik(1094,INDEL::INDEL, 5, insertSeq10);

        const std::pair<int,int> res = test_end_pin_indel_placement(ik);
        BOOST_CHECK_EQUAL(res.first,1004);
        BOOST_CHECK_EQUAL(res.second,0);
    }

    {
        // trailing edge swap:
        IndelKey ik(1095,INDEL::INDEL, 5, insertSeq10);

        const std::pair<int,int> res = test_end_pin_indel_placement(ik);
        BOOST_CHECK_EQUAL(res.first,994);
        BOOST_CHECK_EQUAL(res.second,0);
    }

    {
        // trailing edge swap:
        IndelKey ik(1095,INDEL::INDEL, 5, insertSeq10);

        const std::pair<int,int> res = test_end_pin_indel_placement(ik,95);
        BOOST_CHECK_EQUAL(res.first,999);
        BOOST_CHECK_EQUAL(res.second,0);
    }

    {
        // interfering indel should cause an exception:
        const IndelKey ik(1074,INDEL::INDEL,1);
        BOOST_CHECK_THROW(test_end_pin_indel_placement(ik),blt_exception);
    }
}


BOOST_AUTO_TEST_SUITE_END()
