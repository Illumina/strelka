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

#include "alignment_util.hh"

//#define DEBUG_AU_TEST

#ifdef DEBUG_AU_TEST
#include <iostream>
namespace
{
std::ostream& log_os(std::cerr);
}
#endif



BOOST_AUTO_TEST_SUITE( test_alignment_util )


BOOST_AUTO_TEST_CASE( test_alignment_range_types )
{
    alignment al;
    al.pos = 100;
    ALIGNPATH::cigar_to_apath("2S2I2M2I2S",al.path);

    known_pos_range pr1 = getStrictAlignmentRange(al);
    BOOST_REQUIRE_EQUAL(pr1,known_pos_range(100,102));

    known_pos_range pr2 = get_soft_clip_alignment_range(al);
    BOOST_REQUIRE_EQUAL(pr2,known_pos_range(98,104));

    known_pos_range pr3 = get_alignment_range(al);
    BOOST_REQUIRE_EQUAL(pr3,known_pos_range(96,106));
}

BOOST_AUTO_TEST_CASE( test_alignment_zone )
{
    alignment al;
    al.pos = 100;
    ALIGNPATH::cigar_to_apath("2S2I2M2I2S",al.path);

    known_pos_range pr2 = get_alignment_zone(al,20);
    BOOST_REQUIRE_EQUAL(pr2,known_pos_range(86,116));
}


static
alignment
get_test_alignment()
{
    alignment al;
    al.pos=1000;
    al.is_fwd_strand=true;
    ALIGNPATH::cigar_to_apath("100M",al.path);
    return al;
}


BOOST_AUTO_TEST_CASE( test_remove_edge_deletion )
{

    alignment al(get_test_alignment());

    alignment al2(remove_edge_deletions(al,true,true));

    BOOST_CHECK_EQUAL(al, al2);

    ALIGNPATH::cigar_to_apath("3D100M",al2.path);

    alignment al3(remove_edge_deletions(al2,true,true));

    BOOST_CHECK_EQUAL(al.path, al3.path);

    BOOST_CHECK_EQUAL(static_cast<int>(al3.pos), 1003);
}



BOOST_AUTO_TEST_CASE( test_matchify_insertions )
{
    alignment al;
    al.pos = 100;
    ALIGNPATH::cigar_to_apath("2S2I2M2I2M2I2S",al.path);

    const alignment ali1 = matchify_edge_insertions(al,true,false);
    const alignment ali2 = matchify_edge_insertions(al,false,true);
    const alignment ali12 = matchify_edge_insertions(al,true,true);

    BOOST_REQUIRE_EQUAL(ali1.pos,98);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(ali1.path),"2S4M2I2M2I2S");

    BOOST_REQUIRE_EQUAL(ali2.pos,100);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(ali2.path),"2S2I2M2I4M2S");

    BOOST_REQUIRE_EQUAL(ali12.pos,98);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(ali12.path),"2S4M2I4M2S");
}


BOOST_AUTO_TEST_CASE( test_matchify_soft_clip )
{
    alignment al;
    al.pos = 100;
    ALIGNPATH::cigar_to_apath("2S2I2M2I2M2S",al.path);

    const alignment als = matchify_edge_soft_clip(al);

    BOOST_REQUIRE_EQUAL(als.pos,98);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(als.path),"2M2I2M2I4M");
}


BOOST_AUTO_TEST_CASE( test_getLowestFwdReadPosForRefRange )
{
    alignment al;
    al.pos = 100;
    ALIGNPATH::cigar_to_apath("10M",al.path);

    BOOST_REQUIRE_EQUAL(getLowestFwdReadPosForRefRange(al, known_pos_range(100, 101)),0);

    ALIGNPATH::cigar_to_apath("2S2I2M2I2M2S",al.path);

    BOOST_REQUIRE_EQUAL(getLowestFwdReadPosForRefRange(al, known_pos_range(100, 103)),4);
    BOOST_REQUIRE_EQUAL(getLowestFwdReadPosForRefRange(al, known_pos_range(102, 103)),8);
    BOOST_REQUIRE_EQUAL(getLowestFwdReadPosForRefRange(al, known_pos_range(99, 99)),-1);
    BOOST_REQUIRE_EQUAL(getLowestFwdReadPosForRefRange(al, known_pos_range(110, 110)),-1);

    al.is_fwd_strand = false;

    BOOST_REQUIRE_EQUAL(getLowestFwdReadPosForRefRange(al, known_pos_range(100, 101)),7);
    BOOST_REQUIRE_EQUAL(getLowestFwdReadPosForRefRange(al, known_pos_range(100, 103)),3);
}


BOOST_AUTO_TEST_SUITE_END()
