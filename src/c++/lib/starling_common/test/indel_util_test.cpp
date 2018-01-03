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

#include "indel_util.hh"


BOOST_AUTO_TEST_SUITE( indel_util )


BOOST_AUTO_TEST_CASE( test_indel_intersect_deletion )
{
    // 5-base deletion starting at zero-indexed position 9:
    const IndelKey ikdel1(9,INDEL::INDEL,5);

    // left-side deletion:
    BOOST_CHECK(! is_indel_conflict(ikdel1,IndelKey(7,INDEL::INDEL,1)) );
    BOOST_CHECK(  is_indel_conflict(ikdel1,IndelKey(8,INDEL::INDEL,1)) );

    // right-side deletion:
    BOOST_CHECK(  is_indel_conflict(ikdel1,IndelKey(14,INDEL::INDEL,1)) );
    BOOST_CHECK(! is_indel_conflict(ikdel1,IndelKey(15,INDEL::INDEL,1)) );

    // super-set deletion
    BOOST_CHECK(  is_indel_conflict(ikdel1,IndelKey(4,INDEL::INDEL,15)) );

    // sub-set deletion
    BOOST_CHECK(  is_indel_conflict(ikdel1,IndelKey(12,INDEL::INDEL,1)) );
}



static
bool
test_single_intersect_case(
    const IndelKey& indelKey,
    const bool expect)
{
    // test case represents a 10-base alignment from 10-20 (1-indexed range)
    static const known_pos_range tpr(9,20);

    return (expect==is_range_intersect_indel_breakpoints(tpr,indelKey));
}


BOOST_AUTO_TEST_CASE( test_range_intersect_delete )
{
    static const INDEL::index_t itype(INDEL::INDEL);

    // simple interior case:
    BOOST_CHECK(test_single_intersect_case(IndelKey(14,itype,1),true));

    // left-side edge cases:
    BOOST_CHECK(test_single_intersect_case(IndelKey(8,itype,1),false));
    BOOST_CHECK(test_single_intersect_case(IndelKey(9,itype,1),true));

    // right-side edge cases:
    BOOST_CHECK(test_single_intersect_case(IndelKey(19,itype,1),true));
    BOOST_CHECK(test_single_intersect_case(IndelKey(20,itype,1),false));
}


BOOST_AUTO_TEST_CASE( test_range_intersect_insert )
{
    static const INDEL::index_t itype(INDEL::INDEL);
    static const char* insertSeq = "A";

    // simple interior case:
    BOOST_CHECK(test_single_intersect_case(IndelKey(14,itype, 0, insertSeq),true));

    // left-side edge cases:
    BOOST_CHECK(test_single_intersect_case(IndelKey(9,itype, 0, insertSeq),false));
    BOOST_CHECK(test_single_intersect_case(IndelKey(10,itype, 0, insertSeq),true));

    // right-side edge cases:
    BOOST_CHECK(test_single_intersect_case(IndelKey(19,itype, 0, insertSeq),true));
    BOOST_CHECK(test_single_intersect_case(IndelKey(20,itype, 0, insertSeq),false));
}

BOOST_AUTO_TEST_CASE( test_range_intersect_swap )
{
    static const INDEL::index_t itype(INDEL::INDEL);
    static const char* insertSeq = "A";

    // simple interior case:
    BOOST_CHECK(test_single_intersect_case(IndelKey(14,itype,1,insertSeq),true));

    // left-side edge cases:
    BOOST_CHECK(test_single_intersect_case(IndelKey(8,itype,1, insertSeq),false));
    BOOST_CHECK(test_single_intersect_case(IndelKey(9,itype,1, insertSeq),true));

    // right-side edge cases:
    BOOST_CHECK(test_single_intersect_case(IndelKey(19,itype,1, insertSeq),true));
    BOOST_CHECK(test_single_intersect_case(IndelKey(20,itype,1, insertSeq),false));
}



static
bool
test_single_adjacent_case(
    const IndelKey& indelKey,
    const bool expect)
{
    // test case represents a 10-base alignment from 10-20 (1-indexed range)
    static const known_pos_range tpr(9,20);

    return (expect==is_range_adjacent_indel_breakpoints(tpr,indelKey));
}


BOOST_AUTO_TEST_CASE( test_range_adjacent_delete )
{
    static const INDEL::index_t itype(INDEL::INDEL);

    // simple interior case:
    BOOST_CHECK(test_single_adjacent_case(IndelKey(14,itype,1),true));

    // left-side edge cases:
    BOOST_CHECK(test_single_adjacent_case(IndelKey(7,itype,1),false));
    BOOST_CHECK(test_single_adjacent_case(IndelKey(8,itype,1),true));

    // right-side edge cases:
    BOOST_CHECK(test_single_adjacent_case(IndelKey(20,itype,1),true));
    BOOST_CHECK(test_single_adjacent_case(IndelKey(21,itype,1),false));
}


BOOST_AUTO_TEST_CASE( test_range_adjacent_insert )
{
    static const INDEL::index_t itype(INDEL::INDEL);
    static const char* insertSeq = "A";

    // simple interior case:
    BOOST_CHECK(test_single_adjacent_case(IndelKey(14,itype, 0, insertSeq),true));

    // left-side edge cases:
    BOOST_CHECK(test_single_adjacent_case(IndelKey(8,itype, 0, insertSeq),false));
    BOOST_CHECK(test_single_adjacent_case(IndelKey(9,itype, 0, insertSeq),true));

    // right-side edge cases:
    BOOST_CHECK(test_single_adjacent_case(IndelKey(20,itype, 0, insertSeq),true));
    BOOST_CHECK(test_single_adjacent_case(IndelKey(21,itype, 0, insertSeq),false));
}

BOOST_AUTO_TEST_CASE( test_range_adjacent_swap )
{
    static const INDEL::index_t itype(INDEL::INDEL);
    static const char* insertSeq = "A";

    // simple interior case:
    BOOST_CHECK(test_single_adjacent_case(IndelKey(14,itype,1, insertSeq),true));

    // left-side edge cases:
    BOOST_CHECK(test_single_adjacent_case(IndelKey(7,itype,1, insertSeq),false));
    BOOST_CHECK(test_single_adjacent_case(IndelKey(8,itype,1, insertSeq),true));

    // right-side edge cases:
    BOOST_CHECK(test_single_adjacent_case(IndelKey(20,itype,1, insertSeq),true));
    BOOST_CHECK(test_single_adjacent_case(IndelKey(21,itype,1, insertSeq),false));
}

BOOST_AUTO_TEST_SUITE_END()

