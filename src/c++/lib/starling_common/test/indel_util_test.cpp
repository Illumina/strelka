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

#include "indel_util.hh"


BOOST_AUTO_TEST_SUITE( indel_util )


BOOST_AUTO_TEST_CASE( test_indel_intersect_deletion )
{

    // 5-base deletion starting at zero-indexed position 9:
    const indel_key ikdel1(9,INDEL::DELETE,5);

    // left-side deletion:
    BOOST_CHECK(! is_indel_conflict(ikdel1,indel_key(7,INDEL::DELETE,1)) );
    BOOST_CHECK(  is_indel_conflict(ikdel1,indel_key(8,INDEL::DELETE,1)) );

    // right-side deletion:
    BOOST_CHECK(  is_indel_conflict(ikdel1,indel_key(14,INDEL::DELETE,1)) );
    BOOST_CHECK(! is_indel_conflict(ikdel1,indel_key(15,INDEL::DELETE,1)) );

    // super-set deletion
    BOOST_CHECK(  is_indel_conflict(ikdel1,indel_key(4,INDEL::DELETE,15)) );

    // sub-set deletion
    BOOST_CHECK(  is_indel_conflict(ikdel1,indel_key(12,INDEL::DELETE,1)) );
}



static
bool
test_single_intersect_case(const indel_key& ik,
                           const bool expect)
{

    // test case represents a 10-base alignment from 10-20 (1-indexed range)
    static const known_pos_range tpr(9,20);

    return (expect==is_range_intersect_indel_breakpoints(tpr,ik));
}


BOOST_AUTO_TEST_CASE( test_range_intersect_delete )
{

    static const INDEL::index_t itype(INDEL::DELETE);

    // simple interior case:
    BOOST_CHECK(test_single_intersect_case(indel_key(14,itype,1),true));

    // left-side edge cases:
    BOOST_CHECK(test_single_intersect_case(indel_key(8,itype,1),false));
    BOOST_CHECK(test_single_intersect_case(indel_key(9,itype,1),true));

    // right-side edge cases:
    BOOST_CHECK(test_single_intersect_case(indel_key(19,itype,1),true));
    BOOST_CHECK(test_single_intersect_case(indel_key(20,itype,1),false));
}


BOOST_AUTO_TEST_CASE( test_range_intersect_insert )
{

    static const INDEL::index_t itype(INDEL::INSERT);

    // simple interior case:
    BOOST_CHECK(test_single_intersect_case(indel_key(14,itype,1),true));

    // left-side edge cases:
    BOOST_CHECK(test_single_intersect_case(indel_key(9,itype,1),false));
    BOOST_CHECK(test_single_intersect_case(indel_key(10,itype,1),true));

    // right-side edge cases:
    BOOST_CHECK(test_single_intersect_case(indel_key(19,itype,1),true));
    BOOST_CHECK(test_single_intersect_case(indel_key(20,itype,1),false));
}

BOOST_AUTO_TEST_CASE( test_range_intersect_swap )
{

    static const INDEL::index_t itype(INDEL::SWAP);

    // simple interior case:
    BOOST_CHECK(test_single_intersect_case(indel_key(14,itype,1,1),true));

    // left-side edge cases:
    BOOST_CHECK(test_single_intersect_case(indel_key(8,itype,1,1),false));
    BOOST_CHECK(test_single_intersect_case(indel_key(9,itype,1,1),true));

    // right-side edge cases:
    BOOST_CHECK(test_single_intersect_case(indel_key(19,itype,1,1),true));
    BOOST_CHECK(test_single_intersect_case(indel_key(20,itype,1,1),false));
}



static
bool
test_single_adjacent_case(const indel_key& ik,
                          const bool expect)
{

    // test case represents a 10-base alignment from 10-20 (1-indexed range)
    static const known_pos_range tpr(9,20);

    return (expect==is_range_adjacent_indel_breakpoints(tpr,ik));
}


BOOST_AUTO_TEST_CASE( test_range_adjacent_delete )
{

    static const INDEL::index_t itype(INDEL::DELETE);

    // simple interior case:
    BOOST_CHECK(test_single_adjacent_case(indel_key(14,itype,1),true));

    // left-side edge cases:
    BOOST_CHECK(test_single_adjacent_case(indel_key(7,itype,1),false));
    BOOST_CHECK(test_single_adjacent_case(indel_key(8,itype,1),true));

    // right-side edge cases:
    BOOST_CHECK(test_single_adjacent_case(indel_key(20,itype,1),true));
    BOOST_CHECK(test_single_adjacent_case(indel_key(21,itype,1),false));
}


BOOST_AUTO_TEST_CASE( test_range_adjacent_insert )
{

    static const INDEL::index_t itype(INDEL::INSERT);

    // simple interior case:
    BOOST_CHECK(test_single_adjacent_case(indel_key(14,itype,1),true));

    // left-side edge cases:
    BOOST_CHECK(test_single_adjacent_case(indel_key(8,itype,1),false));
    BOOST_CHECK(test_single_adjacent_case(indel_key(9,itype,1),true));

    // right-side edge cases:
    BOOST_CHECK(test_single_adjacent_case(indel_key(20,itype,1),true));
    BOOST_CHECK(test_single_adjacent_case(indel_key(21,itype,1),false));
}

BOOST_AUTO_TEST_CASE( test_range_adjacent_swap )
{

    static const INDEL::index_t itype(INDEL::SWAP);

    // simple interior case:
    BOOST_CHECK(test_single_adjacent_case(indel_key(14,itype,1,1),true));

    // left-side edge cases:
    BOOST_CHECK(test_single_adjacent_case(indel_key(7,itype,1,1),false));
    BOOST_CHECK(test_single_adjacent_case(indel_key(8,itype,1,1),true));

    // right-side edge cases:
    BOOST_CHECK(test_single_adjacent_case(indel_key(20,itype,1,1),true));
    BOOST_CHECK(test_single_adjacent_case(indel_key(21,itype,1,1),false));
}

BOOST_AUTO_TEST_SUITE_END()

