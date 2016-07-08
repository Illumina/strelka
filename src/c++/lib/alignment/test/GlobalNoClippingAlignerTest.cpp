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

///
/// \author Sangtae Kim
///

#include "boost/test/unit_test.hpp"

#include "GlobalNoClippingAligner.hh"

#include "blt_util/align_path.hh"

#include <string>



BOOST_AUTO_TEST_SUITE( test_GlobalAligner )

typedef short int score_t;

static
AlignmentResult<score_t>
testAlign(
    const std::string& seq,
    const std::string& ref)
{
    AlignmentScores<score_t> scores(1,-4,-5,-1,-1000);
    GlobalNoClippingAligner<score_t> aligner(scores);
    AlignmentResult<score_t> result;
    aligner.align(seq.begin(),seq.end(),ref.begin(),ref.end(),result);

    return result;
}


BOOST_AUTO_TEST_CASE(test_GlobalAlignerLeadingDeletion )
{
    static const std::string seq("CDEF");
    static const std::string ref("ABCDEF");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    // leading deletion is recorded not in apth but in beginPos
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"4=");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,2);
}


BOOST_AUTO_TEST_CASE( test_GlobalAlignerTrailingDeletetion )
{
    static const std::string seq("ABCD");
    static const std::string ref("ABCDEF");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    // trailing deletion is recorded in apath
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"4=2D");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,0);
}

BOOST_AUTO_TEST_CASE(test_GlobalAlignerLeadingInsertion )
{
    static const std::string seq("ABCDEF");
    static const std::string ref("CDEF");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    // leading insertion is recorded as a soft clip
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"2S4=");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,0);
}


BOOST_AUTO_TEST_CASE( test_GlobalAlignerTrailingInsertion )
{
    static const std::string seq("ABCDEF");
    static const std::string ref("ABCD");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    // trailing insertion is recorded in apath
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"4=2I");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,0);
}



BOOST_AUTO_TEST_CASE( test_GlobalAlignerInsertDelete2 )
{
    static const std::string seq("BBBBBBCDEXYHIKLMMMM");
    static const std::string ref("ABBBBBBCDEFGHIKLMMMMN");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"9=2X8=1D");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,1);
}


BOOST_AUTO_TEST_SUITE_END()

