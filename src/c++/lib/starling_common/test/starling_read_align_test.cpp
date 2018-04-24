//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
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

// required to mock-up a read segment:
#include "htsapi/align_path_bam_util.hh"
#include "starling_common/starling_read.hh"
#include "test/starling_base_options_test.hh"



BOOST_AUTO_TEST_SUITE( starling_read_align )

/// \brief Return the alignment resulting from inserting a single indel into a standard background:
static
CandidateAlignment
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
        CandidateAlignment cal = test_indel_placement(ik);
        path_compare("50M10D15M1D35M",cal.al.path);
    }

    {
        // basic insert
        IndelKey ik(1050,INDEL::INDEL, 0, insertSeq10);
        CandidateAlignment cal = test_indel_placement(ik);
        path_compare("50M10I25M1D15M",cal.al.path);
    }

    {
        // basic swap
        IndelKey ik(1050,INDEL::INDEL, 5, insertSeq10);
        CandidateAlignment cal = test_indel_placement(ik);
        path_compare("50M5D10I20M1D20M",cal.al.path);
    }


    {
        // trailing edge insert
        IndelKey ik(1091,INDEL::INDEL, 0, insertSeq10);
        CandidateAlignment cal = test_indel_placement(ik);
        path_compare("75M1D15M10I",cal.al.path);
        BOOST_CHECK_EQUAL(cal.trailing_indel_key,ik);
    }

    {
        // trailing edge insert
        IndelKey ik(1096,INDEL::INDEL, 0, insertSeq10);
        CandidateAlignment cal = test_indel_placement(ik);
        path_compare("75M1D20M5I",cal.al.path);
        BOOST_CHECK_EQUAL(cal.trailing_indel_key,ik);
    }

    {
        // trailing edge insert miss
        IndelKey ik(1101,INDEL::INDEL, 0, insertSeq10);
        CandidateAlignment cal = test_indel_placement(ik);
        path_compare("75M1D25M",cal.al.path);
        BOOST_CHECK_EQUAL(cal.trailing_indel_key.type,INDEL::NONE);
    }

    {
        // leading edge insert miss
        IndelKey ik(1000,INDEL::INDEL, 0, insertSeq10);
        CandidateAlignment cal = test_indel_placement(ik);
        path_compare("75M1D25M",cal.al.path);
        BOOST_CHECK_EQUAL(cal.leading_indel_key.type,INDEL::NONE);
        BOOST_CHECK_EQUAL(cal.al.pos,1000);
    }

    {
        // leading edge insert
        IndelKey ik(1000,INDEL::INDEL, 0, insertSeq10);
        CandidateAlignment cal = test_indel_placement(ik,5);
        path_compare("5I75M1D20M",cal.al.path);
        BOOST_CHECK_EQUAL(cal.leading_indel_key,ik);
        BOOST_CHECK_EQUAL(cal.al.pos,1000);
    }

    {
        // leading edge swap
        IndelKey ik(1000,INDEL::INDEL, 5, insertSeq10);
        CandidateAlignment cal = test_indel_placement(ik,5);
        path_compare("5I5D70M1D25M",cal.al.path);
        BOOST_CHECK_EQUAL(cal.leading_indel_key,ik);
        BOOST_CHECK_EQUAL(cal.al.pos,1000);
    }

    {
        // leading edge delete:
        IndelKey ik(1000,INDEL::INDEL,10);
        CandidateAlignment cal = test_indel_placement(ik);
        path_compare("10D65M1D35M",cal.al.path);
        BOOST_CHECK_EQUAL(cal.leading_indel_key,ik);
        BOOST_CHECK_EQUAL(cal.al.pos,1000);
    }

    {
        // trailing edge delete:
        IndelKey ik(1101,INDEL::INDEL,10);
        CandidateAlignment cal = test_indel_placement(ik);
        path_compare("75M1D25M10D",cal.al.path);
        BOOST_CHECK_EQUAL(cal.trailing_indel_key,ik);
        BOOST_CHECK_EQUAL(cal.al.pos,1000);
    }

    {
        // trailing off-edge delete:
        IndelKey ik(1102,INDEL::INDEL,10);
        CandidateAlignment cal = test_indel_placement(ik);
        path_compare("75M1D25M",cal.al.path);
        BOOST_CHECK_EQUAL(cal.leading_indel_key.type,INDEL::NONE);
        BOOST_CHECK_EQUAL(cal.al.pos,1000);
    }
}



/// \brief Test results of placing a single indel into a standard background.
///
/// \return A pair of positions corresponding to: (ref_start_pos,read_start_pos)
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



BOOST_AUTO_TEST_CASE( test_realign_and_score_read )
{
    // the realigner normally throws away edge soft-clipping to try to produce more indel information, this
    // is a test to ensure that in situations where it makes no sense to force an edge-softclip to the match
    // state, that the realigner does not do so (see STREL-459)
    //
    // this also demos the full procedure required to mock realignAndScoreRead, which suggest a number of
    // helpful test utility functions -- in particular generating a bam_record object directly from a valid
    // SAM string
    {
        starling_base_options_test opt;
        opt.isRetainOptimalSoftClipping = true;

        starling_base_deriv_options dopt(opt);
        starling_sample_options sample_opt(opt);
        reference_contig_segment ref;
        ref.seq() = "ACGTACGTACGTACGTACGT";

        known_pos_range realign_buffer_range(0, 20);

        const unsigned sampleIndex(0);

        // to mock up the indel buffer with one candidate indel, we have to initialize it with one fake sample,
        // finalize, and then insert an indel observation
        IndelBuffer indelBuffer(opt, dopt, ref);
        depth_buffer db;
        depth_buffer db2;
        indelBuffer.registerSample(db, db2, false);
        indelBuffer.finalizeSamples();
        {
            IndelObservation obs;
            obs.key = IndelKey(4, INDEL::INDEL, 1);
            obs.data.is_external_candidate = true;

            indelBuffer.addIndelObservation(sampleIndex, obs);
        }

        // the read segment is the hardest piece of this to mock up:
        // 1) mock up the underlying bam record:
        bam_record bamRead;
        bamRead.set_qname("FOOREAD");
        const char read[] = "GTACGG";
        const uint8_t qual[] = {40, 40, 40, 40, 40, 40};
        bamRead.set_readqual(read, qual);

        // set alignment
        alignment al;
        al.pos = 2;
        ALIGNPATH::cigar_to_apath("5M1S", al.path);

        // set bam record from alignment
        bam1_t& br(*(bamRead.get_data()));
        br.core.pos = al.pos;
        edit_bam_cigar(al.path, br);

        // 2) mock up the starling read
        starling_read sread(bamRead, al, MAPLEVEL::UNKNOWN, 0);

        // 3) finally, get read_segment from starling_read
        read_segment& rseg(sread.get_full_segment());

        realignAndScoreRead(opt, dopt, sample_opt, ref, realign_buffer_range, sampleIndex, rseg,
                            indelBuffer);

        BOOST_REQUIRE(not rseg.is_realigned);
    }
}

BOOST_AUTO_TEST_SUITE_END()
