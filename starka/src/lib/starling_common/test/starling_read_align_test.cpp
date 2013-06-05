// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

#include "boost/test/unit_test.hpp"

#include "starling_read_align.cpp"



BOOST_AUTO_TEST_SUITE( starling_read_align )

// test results of placing a single indel into a standard background:
static
candidate_alignment
test_indel_placement(const indel_key ik,
                     const pos_t read_start_pos = 0) {

    const pos_t ref_start_pos(1000);
    const bool is_fwd_strand(true);
    const unsigned read_length(100);

    indel_set_t iset;
    const indel_key ikfixed(1075,INDEL::DELETE,1);
    iset.insert(ikfixed);
    iset.insert(ik);

    return make_start_pos_alignment(ref_start_pos,read_start_pos,is_fwd_strand,read_length,iset);
}


static
void
path_compare(const std::string& expect,
             const ALIGNPATH::path_t& result) {

    using namespace ALIGNPATH;

    path_t tmp_path;
    cigar_to_apath(expect.c_str(),tmp_path);
    BOOST_CHECK_EQUAL(tmp_path,result);
}


BOOST_AUTO_TEST_CASE( test_make_start_pos_alignment ) {

    { // basic delete
        indel_key ik(1050,INDEL::DELETE,10);
        candidate_alignment cal = test_indel_placement(ik);
        path_compare("50M10D15M1D35M",cal.al.path);
    }

    { // basic insert
        indel_key ik(1050,INDEL::INSERT,10);
        candidate_alignment cal = test_indel_placement(ik);
        path_compare("50M10I25M1D15M",cal.al.path);
    }

    { // basic swap
        indel_key ik(1050,INDEL::SWAP,10,5);
        candidate_alignment cal = test_indel_placement(ik);
        path_compare("50M5D10I20M1D20M",cal.al.path);
    }


    { // trailing edge insert
        indel_key ik(1091,INDEL::INSERT,10);
        candidate_alignment cal = test_indel_placement(ik);
        path_compare("75M1D15M10I",cal.al.path);
        BOOST_CHECK_EQUAL(cal.trailing_indel_key,ik);
    }

    { // trailing edge insert
        indel_key ik(1096,INDEL::INSERT,10);
        candidate_alignment cal = test_indel_placement(ik);
        path_compare("75M1D20M5I",cal.al.path);
        BOOST_CHECK_EQUAL(cal.trailing_indel_key,ik);
    }

    { // trailing edge insert miss
        indel_key ik(1101,INDEL::INSERT,10);
        candidate_alignment cal = test_indel_placement(ik);
        path_compare("75M1D25M",cal.al.path);
        BOOST_CHECK_EQUAL(cal.trailing_indel_key.type,INDEL::NONE);
    }

    { // leading edge insert miss
        indel_key ik(1000,INDEL::INSERT,10);
        candidate_alignment cal = test_indel_placement(ik);
        path_compare("75M1D25M",cal.al.path);
        BOOST_CHECK_EQUAL(cal.leading_indel_key.type,INDEL::NONE);
        BOOST_CHECK_EQUAL(cal.al.pos,1000);
    }

    { // leading edge insert
        indel_key ik(1000,INDEL::INSERT,10);
        candidate_alignment cal = test_indel_placement(ik,5);
        path_compare("5I75M1D20M",cal.al.path);
        BOOST_CHECK_EQUAL(cal.leading_indel_key,ik);
        BOOST_CHECK_EQUAL(cal.al.pos,1000);
    }

    { // leading edge swap
        indel_key ik(1000,INDEL::SWAP,10,5);
        candidate_alignment cal = test_indel_placement(ik,5);
        path_compare("5I5D70M1D25M",cal.al.path);
        BOOST_CHECK_EQUAL(cal.leading_indel_key,ik);
        BOOST_CHECK_EQUAL(cal.al.pos,1000);
    }
}



BOOST_AUTO_TEST_CASE( test_end_pin_start_pos ) {

    const pos_t ref_end_pos(1100);
    const pos_t read_end_pos(100);
    const unsigned read_length(100);

    indel_key ik(1050,INDEL::DELETE,10);
    indel_set_t iset;
    iset.insert(ik);

    pos_t new_ref_start_pos(0);
    pos_t new_read_start_pos(0);

    get_end_pin_start_pos(iset,read_length,ref_end_pos,read_end_pos,new_ref_start_pos,new_read_start_pos);

    BOOST_CHECK_EQUAL(static_cast<int>(new_ref_start_pos),990);
    BOOST_CHECK_EQUAL(static_cast<int>(new_read_start_pos),0);
}

BOOST_AUTO_TEST_SUITE_END()
