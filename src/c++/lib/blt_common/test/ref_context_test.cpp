// -*- mode: c++; indent-tabs-mode: nil; -*-
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

#include "boost/test/unit_test.hpp"

#include "ref_context.hh"


BOOST_AUTO_TEST_SUITE( ref_context )


BOOST_AUTO_TEST_CASE( test_left_shifted_hpol_size )
{
    reference_contig_segment ref;
    ref.seq() = "TTGTTTGAGATTT";

    BOOST_REQUIRE_EQUAL(get_left_shifted_hpol_size(0,ref),2u);
    BOOST_REQUIRE_EQUAL(get_left_shifted_hpol_size(1,ref),1u);
    BOOST_REQUIRE_EQUAL(get_left_shifted_hpol_size(2,ref),1u);
    BOOST_REQUIRE_EQUAL(get_left_shifted_hpol_size(3,ref),3u);
    BOOST_REQUIRE_EQUAL(get_left_shifted_hpol_size(4,ref),1u);
}


static
void
single_snp_hpol_test(const pos_t pos,
                     const reference_contig_segment& ref,
                     const unsigned expect)
{

    const unsigned result(get_snp_hpol_size(pos,ref));
    BOOST_CHECK_EQUAL(result,expect);
}


BOOST_AUTO_TEST_CASE( test_snp_hpol_size )
{
    reference_contig_segment ref;
    ref.seq() = "TGTTTGAGATTT";

    single_snp_hpol_test(0,ref,2);
    single_snp_hpol_test(1,ref,5);
}



static
void
single_ihpol_test(const pos_t pos,
                  const reference_contig_segment& ref,
                  const unsigned expect)
{

    const unsigned result(get_interrupted_hpol_size(pos,ref));
    BOOST_CHECK_EQUAL(result,expect);
}


BOOST_AUTO_TEST_CASE( test_interrupted_hpol_size )
{
    reference_contig_segment ref;
    ref.seq() = "TGTTTGAGATTT";

    single_ihpol_test(0,ref,4);
    single_ihpol_test(1,ref,4);
}


BOOST_AUTO_TEST_SUITE_END()

