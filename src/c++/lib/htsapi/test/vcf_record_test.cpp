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

#include "htsapi/vcf_record.hh"
#include "blt_util/reference_contig_segment.hh"

#include "boost/test/unit_test.hpp"


BOOST_AUTO_TEST_SUITE( vcf_record_test )

BOOST_AUTO_TEST_CASE( test_snv )
{
    reference_contig_segment ref_contig;
    ref_contig.seq().assign("AGTCCAC");

    vcf_record vcfr;
    vcfr.set("chr1\t1\t.\tA\tT\n", ref_contig);
    BOOST_REQUIRE(not vcfr.is_indel());
    BOOST_REQUIRE(vcfr.is_snv());
    BOOST_REQUIRE(not vcfr.is_ref_site());
}

BOOST_AUTO_TEST_CASE( test_site )
{
    reference_contig_segment ref_contig;
    ref_contig.seq().assign("AGTCCAC");

    vcf_record vcfr;
    vcfr.set("chr1\t1\t.\tA\t.\n", ref_contig);
    BOOST_REQUIRE(not vcfr.is_indel());
    BOOST_REQUIRE(not vcfr.is_snv());
    BOOST_REQUIRE(vcfr.is_ref_site());
}


BOOST_AUTO_TEST_SUITE_END()

