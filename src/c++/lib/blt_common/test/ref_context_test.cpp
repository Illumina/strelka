// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

#include "boost/test/unit_test.hpp"

#include "ref_context.hh"


BOOST_AUTO_TEST_SUITE( ref_context )


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

