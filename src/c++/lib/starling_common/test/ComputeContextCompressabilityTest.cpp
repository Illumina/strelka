//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
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
/// \author Konrad Scheffler
///

#include <blt_util/reference_contig_segment.hh>
#include <starling_common/AlleleReportInfoUtil.hh>
#include "boost/test/unit_test.hpp"

BOOST_AUTO_TEST_SUITE( test_computeContextCompressability )

// Checks whether shortest unencountered sequence is determined correctly
BOOST_AUTO_TEST_CASE( test_shortestUnencountered )
{
    reference_contig_segment ref;
    ref.seq() = "AAACGACGAC";

    // With the 2nd base as the first to be encoded, look to the right. After encountering
    // a given substring, what is the shortest unencoutered string in what follows?
    const pos_t secondpos(1);

    // Encountered nothing, shortest unencountered string in "AAC..." is A:
    BOOST_REQUIRE_EQUAL(shortestUnencountered(ref,secondpos,0), 1);

    // Encountered A, shortest unencountered string in "ACG..." is AC:
    BOOST_REQUIRE_EQUAL(shortestUnencountered(ref,secondpos,1), 2);

    // Encountered AA, shortest unencountered string in "CGA..." is C:
    BOOST_REQUIRE_EQUAL(shortestUnencountered(ref,secondpos,2), 1);

    // Encountered AAC, shortest unencountered string in "GAC..." is G:
    BOOST_REQUIRE_EQUAL(shortestUnencountered(ref,secondpos,3), 1);

    // Encountered AACG, shortest unencountered string in "ACGAC" is ACGA:
    BOOST_REQUIRE_EQUAL(shortestUnencountered(ref,secondpos,4), 4);

    // Encountered AACGA, shortest unencountered string in "CGAC" is CGAC:
    BOOST_REQUIRE_EQUAL(shortestUnencountered(ref,secondpos,5), 4);

    // Encountered AACGAC, every substring of "GAC" has been encountered.
    // Shortest unencountered will be one beyond the end of the string:
    BOOST_REQUIRE_EQUAL(shortestUnencountered(ref,secondpos,6), 4);

    // With the 2nd last base as the first to be encoded, look to the left:
    const pos_t secondlastpos(8);

    // Encountered nothing, shortest unencountered string in "...A" is A:
    BOOST_REQUIRE_EQUAL(shortestUnencountered(ref,secondlastpos,0,true), 1);

    // Encountered A, shortest unencountered string in "...G" is G:
    BOOST_REQUIRE_EQUAL(shortestUnencountered(ref,secondlastpos,1,true), 1);

    // Encountered GA, shortest unencountered string in "...C" is C:
    BOOST_REQUIRE_EQUAL(shortestUnencountered(ref,secondlastpos,2,true), 1);

    // Encountered CGA, shortest unencountered string in "...ACGA" is ACGA:
    BOOST_REQUIRE_EQUAL(shortestUnencountered(ref,secondlastpos,3,true), 4);

    // Encountered ACGA, shortest unencountered string in "...AACG" is AACG:
    BOOST_REQUIRE_EQUAL(shortestUnencountered(ref,secondlastpos,4,true), 4);
}

// Checks whether contextCompressability is computed correctly
BOOST_AUTO_TEST_CASE( test_computeContextCompressability )
{
    reference_contig_segment ref;
    ref.seq() = "AAAAAAAATGC";
    const pos_t pos(8);

    // For an insertion after the last A (left context: homopolymer, right context: TGC),
    // three keywords will encode 1+2+4=7 bases on the left (and 1+1+1=3 bases on the right):
    BOOST_REQUIRE_EQUAL(computeContextCompressability(ref,pos,pos,3), 7);
}

BOOST_AUTO_TEST_SUITE_END()
