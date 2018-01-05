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

#include "GlobalAligner.hh"

#include "blt_util/align_path.hh"

#include <string>



BOOST_AUTO_TEST_SUITE( test_GlobalAligner )

typedef short int score_t;

static
AlignmentResult<score_t>
testAlign(
    const std::string& seq,
    const std::string& ref,
    const int offEdgeScore = -4,
    const int insertDeleteScore = 0,
    const bool isAllowEdgeInsertion = false,
    const bool isRequireEdgeDeletion = false)
{
    AlignmentScores<score_t> scores(2, -4, -5, -1, offEdgeScore, insertDeleteScore, isAllowEdgeInsertion, isRequireEdgeDeletion);
    GlobalAligner<score_t> aligner(scores);
    AlignmentResult<score_t> result;
    aligner.align(seq.begin(),seq.end(),ref.begin(),ref.end(),result);

    return result;
}


BOOST_AUTO_TEST_CASE( test_GlobalAligner1 )
{
    static const std::string seq("D");
    static const std::string ref("ABCDEF");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"1=");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,3);
}



BOOST_AUTO_TEST_CASE( test_GlobalAlignerDelete )
{
    static const std::string seq("BCDEFHIKLM");
    static const std::string ref("ABCDEFGHIKLMN");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"5=1D5=");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,1);
}


BOOST_AUTO_TEST_CASE( test_GlobalAlignerInsert )
{
    static const std::string seq("BCDEFGXHIKLM");
    static const std::string ref("ABCDEFGHIKLMN");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"6=1I5=");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,1);
}


BOOST_AUTO_TEST_CASE( test_GlobalAlignerInsertDelete )
{
    static const std::string seq("BBBBBBCDXYZHIKLMMMM");
    static const std::string ref("ABBBBBBCDEFGHIKLMMMMN");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"8=3I3D8=");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,1);
}

BOOST_AUTO_TEST_CASE( test_GlobalAlignerInsertDelete2 )
{
    static const std::string seq("BBBBBBCDEXYHIKLMMMM");
    static const std::string ref("ABBBBBBCDEFGHIKLMMMMN");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"9=2X8=");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,1);
}


BOOST_AUTO_TEST_CASE( test_GlobalAlignerShortRef1 )
{
    static const std::string seq("ABCD");
    static const std::string ref("BCD");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"1S3=");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,0);
    BOOST_REQUIRE_EQUAL(result.score,2);
}


BOOST_AUTO_TEST_CASE( test_GlobalAlignerShortRef2 )
{
    static const std::string seq("ABCD");
    static const std::string ref("ABC");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"3=1S");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,0);
    BOOST_REQUIRE_EQUAL(result.score,2);
}


BOOST_AUTO_TEST_CASE( test_GlobalAlignerShortRef3 )
{
    static const std::string seq("ABCD");
    static const std::string ref("B");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"1S1=2S");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,0);
    BOOST_REQUIRE_EQUAL(result.score,-10);
}


// show that the method left aligns a deletion within a repeat
//
BOOST_AUTO_TEST_CASE( test_GlobalAlignerLeftShift )
{
    static const std::string seq("ABCDEFFFFFGHIJKL");
    static const std::string ref("ABCDEFFFFFFGHIJKL");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"5=1D11=");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,0);
}

// show that the method left aligns an insertion within a repeat
//
BOOST_AUTO_TEST_CASE( test_GlobalAlignerLeftShift2 )
{
    static const std::string seq("ABCDEFFFFFFFGHIJKL");
    static const std::string ref("ABCDEFFFFFFGHIJKL");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"5=1I12=");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,0);
}


// show that alignment is global over the query
//
BOOST_AUTO_TEST_CASE( test_GlobalAlignerGlocal )
{
    static const std::string seq("AABCC");
    static const std::string ref("ZZBYY");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"2X1=2X");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,0);
}


// show that soft-clipping weights and edge indel blocking are working correctly by setting
// off-edge score to a prohibitively high value.
//
BOOST_AUTO_TEST_CASE( test_GlobalAlignerNoClip )
{
    static const std::string seq("12ABCDEFGHIJ12");
    static const std::string ref(  "ABCDEFGHIJ");

    AlignmentResult<score_t> result = testAlign(seq,ref,-100);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"1X2I8=2I1X");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,0);
}


// show that soft-clipping weights and edge indel blocking are working correctly by setting
// off-edge score to a prohibitively high value, and allowing edge Indels
//
BOOST_AUTO_TEST_CASE( test_GlobalAlignerNoClipEdgeIndel )
{
    static const std::string seq("12ABCDEFGHIJ12");
    static const std::string ref(  "ABCDEFGHIJ");

    AlignmentResult<score_t> result = testAlign(seq,ref,-100, 0, true);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"2I10=2I");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,0);
    BOOST_REQUIRE_EQUAL(result.score,6);
}


//
BOOST_AUTO_TEST_CASE( test_GlobalAlignerNoClipEdgeIndelMini1 )
{
    static const std::string seq("AB");
    static const std::string ref("A");

    AlignmentResult<score_t> result = testAlign(seq,ref,-100, 0, true);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"1=1I");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,0);
    BOOST_REQUIRE_EQUAL(result.score,-4);
}


//
BOOST_AUTO_TEST_CASE( test_GlobalAlignerNoClipEdgeIndelMini2 )
{
    static const std::string seq("AB");
    static const std::string ref( "B");

    AlignmentResult<score_t> result = testAlign(seq,ref,-100, 0, true);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"1I1=");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,0);
    BOOST_REQUIRE_EQUAL(result.score,-4);
}


//
BOOST_AUTO_TEST_CASE( test_GlobalAlignerEdgeDeletionOff )
{
    static const std::string seq(  "CDEFFFGHIJ");
    static const std::string ref("ABCDEFFFGHIJKL");

    AlignmentResult<score_t> result = testAlign(seq,ref,-100, 0, true);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"10=");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,2);
    BOOST_REQUIRE_EQUAL(result.score,20);
}


//
BOOST_AUTO_TEST_CASE( test_GlobalAlignerEdgeDeletionOn )
{
    static const std::string seq(  "CDEFFFGHIJ");
    static const std::string ref("ABCDEFFFGHIJKL");

    AlignmentResult<score_t> result = testAlign(seq,ref,-100, 0, true, true);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"2D10=2D");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,0);
    BOOST_REQUIRE_EQUAL(result.score,6);
}


//
BOOST_AUTO_TEST_CASE( test_GlobalAlignerEdgeDeletionOnMini1 )
{
    static const std::string seq("A");
    static const std::string ref("AC");

    AlignmentResult<score_t> result = testAlign(seq,ref,-100, 0, true, true);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"1=1D");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,0);
    BOOST_REQUIRE_EQUAL(result.score,-4);
}


//
BOOST_AUTO_TEST_CASE( test_GlobalAlignerEdgeDeletionOnMini2 )
{
    static const std::string seq( "C");
    static const std::string ref("AC");

    AlignmentResult<score_t> result = testAlign(seq,ref,-100, 0, true, true);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"1D1=");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,0);
    BOOST_REQUIRE_EQUAL(result.score,-4);
}


// test insert left shifting when an edge is involved
BOOST_AUTO_TEST_CASE( test_GlobalAlignerNoClipEdgeIndelLeftShiftInsert )
{
    static const std::string seq("ATTT");
    static const std::string ref("AT");

    AlignmentResult<score_t> result = testAlign(seq,ref,-100, 0, true, true);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"1=2I1=");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,0);
    BOOST_REQUIRE_EQUAL(result.score,-3);
}


// test delete left shifting when an edge is involved
BOOST_AUTO_TEST_CASE( test_GlobalAlignerNoClipEdgeIndelLeftShiftDelete )
{
    static const std::string seq("AT");
    static const std::string ref("ATTT");

    AlignmentResult<score_t> result = testAlign(seq,ref,-100, 0, true, true);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"1=2D1=");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,0);
    BOOST_REQUIRE_EQUAL(result.score,-3);
}


//
BOOST_AUTO_TEST_CASE( test_GlobalAlignerNoClipEdgeIndelInsertDelete )
{
    static const std::string seq("GCG");
    static const std::string ref("GCCC");

    AlignmentResult<score_t> result = testAlign(seq,ref,-100, -10, true, true);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"1=1D1=1X");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,0);
}


BOOST_AUTO_TEST_SUITE_END()
