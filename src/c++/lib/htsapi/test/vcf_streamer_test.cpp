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

#include "test_config.h"
#include <iostream>

#include "htsapi/vcf_streamer.hh"

#include "common/Exceptions.hh"

#include "boost/test/unit_test.hpp"


BOOST_AUTO_TEST_SUITE( test_vcf_streamer )

static
const char*
getTestVcfPath()
{
    static const std::string testPath(std::string(TEST_DATA_PATH) + "/vcf_streamer_test.vcf.gz");
    return testPath.c_str();
}

const char*
getTestFastaPath()
{
    static const std::string testPath(std::string(TEST_DATA_PATH) + "/vcf_streamer_test.fa");
    return testPath.c_str();
}


BOOST_AUTO_TEST_CASE( test_vcf_streamer_region )
{
    vcf_streamer vcfs(getTestVcfPath(),"chrP:1-100", getTestFastaPath());

    const vcf_record* vptr(nullptr);

    // first check that valid record does not throw an exception
    BOOST_REQUIRE_NO_THROW( vcfs.next() );
    vptr = vcfs.get_record_ptr();
    assert(vptr != nullptr);

    // VCF record test for
    // chrP    1    MantaINS:29:0:0:0:0:0    T    TTTTCTTTC    ...
    // testing variant assignment, normalize and validity checks, position, and ref
    BOOST_REQUIRE(vptr->isSimpleVariantLocus() );
    BOOST_REQUIRE( vptr->is_indel() );
    BOOST_REQUIRE( vptr->is_normalized());
    BOOST_REQUIRE( vptr->is_match_reference());
    BOOST_REQUIRE( ! vptr->is_snv() );

    BOOST_REQUIRE_EQUAL(vptr->pos, 1);
    BOOST_REQUIRE_EQUAL(vptr->ref,"T");
    BOOST_REQUIRE_EQUAL(vptr->alt.size(),1u);
    BOOST_REQUIRE_EQUAL(vptr->alt[0],"TTTTCTTTC");

    BOOST_REQUIRE( vcfs.next() );
    vptr = vcfs.get_record_ptr();
    assert(vptr != nullptr);

    // VCF record test for
    // chrP    2    MantaDEL:41:0:0:0:3:0    CGAATGGAATG    C    ...
    // testing variant assignment, normalize and validity checks, position, and ref
    BOOST_REQUIRE(vptr->isSimpleVariantLocus() );
    BOOST_REQUIRE( vptr->is_indel() );
    BOOST_REQUIRE( vptr->is_normalized());
    BOOST_REQUIRE( vptr->is_match_reference());
    BOOST_REQUIRE( ! vptr->is_snv() );

    BOOST_REQUIRE_EQUAL(vptr->pos, 2);
    BOOST_REQUIRE_EQUAL(vptr->ref,"CGAATGGAATG");
    BOOST_REQUIRE_EQUAL(vptr->alt.size(),1u);
    BOOST_REQUIRE_EQUAL(vptr->alt[0],"C");

    BOOST_REQUIRE( vcfs.next() );
    vptr = vcfs.get_record_ptr();
    assert(vptr != nullptr);

    // VCF record test for
    // chrP    13  MantaDEL:44:0:0:0:1:0   CCCTGGCCAGCAGATCCACCCTGTCTATACTACCTG    C       ...
    // testing variant assignment, normalize and validity checks, position, and ref
    BOOST_REQUIRE(vptr->isSimpleVariantLocus() );
    BOOST_REQUIRE( vptr->is_indel() );
    BOOST_REQUIRE( vptr->is_normalized());
    BOOST_REQUIRE( vptr->is_match_reference());
    BOOST_REQUIRE( ! vptr->is_snv() );

    BOOST_REQUIRE_EQUAL(vptr->pos, 13);
    BOOST_REQUIRE_EQUAL(vptr->ref,"CCCTGGCCAGCAGATCCACCCTGTCTATACTACCTG");
    BOOST_REQUIRE_EQUAL(vptr->alt.size(),1u);
    BOOST_REQUIRE_EQUAL(vptr->alt[0],"C");

    // also check that a valid record returns true
    BOOST_REQUIRE( vcfs.next() );
    vptr = vcfs.get_record_ptr();
    assert(vptr != nullptr);

    // VCF record test for
    // chrP    39  FAKED   A       T,C     ...
    // testing variant assignment, normalize and validity checks, position, alt size, and alt
    BOOST_REQUIRE(vptr->isSimpleVariantLocus() );
    BOOST_REQUIRE( ! vptr->is_indel() );
    BOOST_REQUIRE( vptr->is_snv() );
    BOOST_REQUIRE( vptr->is_normalized());
    BOOST_REQUIRE( vptr->is_match_reference());

    BOOST_REQUIRE_EQUAL(vptr->pos, 39);
    BOOST_REQUIRE_EQUAL(vptr->ref,"A");
    BOOST_REQUIRE_EQUAL(vptr->alt.size(),2u);
    BOOST_REQUIRE_EQUAL(vptr->alt[0],"T");
    BOOST_REQUIRE_EQUAL(vptr->alt[1],"C");

    BOOST_REQUIRE( vcfs.next() );
    vptr = vcfs.get_record_ptr();
    assert(vptr != nullptr);

    // VCF record test for
    // chrP    40  MantaINS:53:0:0:0:1:0   T       TGCCCTTTGGCAGAGCAGGTGTGCTGTGCTG ...
    // testing variant assignment, normalize and validity checks, position, alt size, and alt
    BOOST_REQUIRE(vptr->isSimpleVariantLocus() );
    BOOST_REQUIRE( vptr->is_indel() );
    BOOST_REQUIRE( ! vptr->is_snv() );
    BOOST_REQUIRE( vptr->is_normalized());
    BOOST_REQUIRE( vptr->is_match_reference());

    BOOST_REQUIRE_EQUAL(vptr->pos, 40);
    BOOST_REQUIRE_EQUAL(vptr->ref,"T");
    BOOST_REQUIRE_EQUAL(vptr->alt.size(),1u);
    BOOST_REQUIRE_EQUAL(vptr->alt[0],"TGCCCTTTGGCAGAGCAGGTGTGCTGTGCTG");

    vcfs.resetRegion("chrQ:1-100");

    // tests that we can reset VCF regions
    BOOST_REQUIRE( vcfs.next() );
    vptr = vcfs.get_record_ptr();
    assert(vptr != nullptr);

    // VCF record test for
    // chrQ   1    NORMALIZED    TA    T    ...
    // testing that deletion is reported as normalized
    BOOST_REQUIRE(vptr->isSimpleVariantLocus() );
    BOOST_REQUIRE( vptr->is_indel() );
    BOOST_REQUIRE( ! vptr->is_snv() );
    BOOST_REQUIRE( vptr->is_normalized());
    BOOST_REQUIRE( vptr->is_match_reference());
    BOOST_REQUIRE_EQUAL(vptr->pos, 1);
    BOOST_REQUIRE_EQUAL(vptr->alt.size(),1u);
    BOOST_REQUIRE_EQUAL(vptr->alt[0],"T");

    BOOST_REQUIRE_THROW( vcfs.next(), illumina::common::LogicException );
    vptr = vcfs.get_record_ptr();
    assert(vptr != nullptr);

    // VCF record test for
    // chrQ   3    UNNORMALIZED_DEL    AA    A    ...
    // testing that deletion is not reported as left-shifted
    BOOST_REQUIRE(vptr->isSimpleVariantLocus() );
    BOOST_REQUIRE( vptr->is_indel() );
    BOOST_REQUIRE( ! vptr->is_snv() );
    BOOST_REQUIRE( ! vptr->is_normalized());
    BOOST_REQUIRE( vptr->is_match_reference());
    BOOST_REQUIRE_EQUAL(vptr->pos, 3);
    BOOST_REQUIRE_EQUAL(vptr->alt.size(),1u);
    BOOST_REQUIRE_EQUAL(vptr->alt[0],"A");

    BOOST_REQUIRE_THROW( vcfs.next(), illumina::common::LogicException );
    vptr = vcfs.get_record_ptr();
    assert(vptr != nullptr);

    // VCF record test for
    // chrQ   5    COSM5012    A    AA    ...
    // testing that insertion is not reported as left-shifted
    BOOST_REQUIRE(vptr->isSimpleVariantLocus() );
    BOOST_REQUIRE( vptr->is_indel() );
    BOOST_REQUIRE( ! vptr->is_snv() );
    BOOST_REQUIRE( ! vptr->is_normalized());
    BOOST_REQUIRE( vptr->is_match_reference());
    BOOST_REQUIRE_EQUAL(vptr->pos, 5);
    BOOST_REQUIRE_EQUAL(vptr->alt.size(),1u);
    BOOST_REQUIRE_EQUAL(vptr->alt[0],"AA");

    BOOST_REQUIRE( vcfs.next() );
    vptr = vcfs.get_record_ptr();
    assert(vptr != nullptr);

    // VCF record test for
    // chrQ   6    VALID_MNV_1    TCG    AGT    ...
    // testing that MNV is reported as normalized
    BOOST_REQUIRE(vptr->isSimpleVariantLocus() );
    BOOST_REQUIRE( vptr->is_indel() );
    BOOST_REQUIRE( ! vptr->is_snv() );
    BOOST_REQUIRE( vptr->is_normalized());
    BOOST_REQUIRE( vptr->is_match_reference());
    BOOST_REQUIRE_EQUAL(vptr->pos, 6);
    BOOST_REQUIRE_EQUAL(vptr->alt.size(),1u);
    BOOST_REQUIRE_EQUAL(vptr->alt[0],"AGT");

    BOOST_REQUIRE_THROW( vcfs.next(), illumina::common::LogicException );
    vptr = vcfs.get_record_ptr();
    assert(vptr != nullptr);

    // VCF record test for
    // chrQ   10    INVALID_MNV_1    TCG    AGG    ...
    // testing that MNV is not reported as normalized
    BOOST_REQUIRE(vptr->isSimpleVariantLocus() );
    BOOST_REQUIRE( vptr->is_indel() );
    BOOST_REQUIRE( ! vptr->is_snv() );
    BOOST_REQUIRE( ! vptr->is_normalized());
    BOOST_REQUIRE( vptr->is_match_reference());
    BOOST_REQUIRE_EQUAL(vptr->pos, 10);
    BOOST_REQUIRE_EQUAL(vptr->alt.size(),1u);
    BOOST_REQUIRE_EQUAL(vptr->alt[0],"AGG");

    BOOST_REQUIRE( vcfs.next() );
    vptr = vcfs.get_record_ptr();
    assert(vptr != nullptr);

    // VCF record test for
    // chrQ   13    INVALID_MNV_2   ACG    AGT    ...
    // testing that MNV is reported as normalized, despite first base
    // matching (Manta returns MNV candidates with reference-padding, so
    // this is to maintain consistency with its inputs)
    BOOST_REQUIRE(vptr->isSimpleVariantLocus() );
    BOOST_REQUIRE( vptr->is_indel() );
    BOOST_REQUIRE( ! vptr->is_snv() );
    BOOST_REQUIRE( vptr->is_normalized());
    BOOST_REQUIRE( vptr->is_match_reference());
    BOOST_REQUIRE_EQUAL(vptr->pos, 13);
    BOOST_REQUIRE_EQUAL(vptr->alt.size(),1u);
    BOOST_REQUIRE_EQUAL(vptr->alt[0],"AGT");

    BOOST_REQUIRE_THROW( vcfs.next(), illumina::common::LogicException );
    vptr = vcfs.get_record_ptr();
    assert(vptr != nullptr);

    // VCF record test for
    // chrQ   16    INVALID_SNV   A    A    ...
    // testing that SNV is not reported as normalized
    BOOST_REQUIRE(vptr->isSimpleVariantLocus() );
    BOOST_REQUIRE( ! vptr->is_indel() );
    BOOST_REQUIRE( vptr->is_snv() );
    BOOST_REQUIRE( ! vptr->is_normalized());
    BOOST_REQUIRE( vptr->is_match_reference());
    BOOST_REQUIRE_EQUAL(vptr->pos, 16);
    BOOST_REQUIRE_EQUAL(vptr->alt.size(),1u);
    BOOST_REQUIRE_EQUAL(vptr->alt[0],"A");

    BOOST_REQUIRE_THROW( vcfs.next(), illumina::common::LogicException );
    vptr = vcfs.get_record_ptr();
    assert(vptr != nullptr);

    // VCF record test for
    // chrQ   17    RIGHT_PAD_INDEL    ACCC    AC    ...
    // testing that indel is not reported as normalized
    BOOST_REQUIRE(vptr->isSimpleVariantLocus() );
    BOOST_REQUIRE( vptr->is_indel() );
    BOOST_REQUIRE( ! vptr->is_snv() );
    BOOST_REQUIRE( ! vptr->is_normalized());
    BOOST_REQUIRE( vptr->is_match_reference());
    BOOST_REQUIRE_EQUAL(vptr->pos, 17);
    BOOST_REQUIRE_EQUAL(vptr->alt.size(),1u);
    BOOST_REQUIRE_EQUAL(vptr->alt[0],"AC");

    BOOST_REQUIRE( vcfs.next() );
    vptr = vcfs.get_record_ptr();
    assert(vptr != nullptr);

    // VCF record test for
    // chrQ   21    COMPLEX    GAGCTGTG   AGCT    ...
    // testing that complex allele is reported as normalized
    BOOST_REQUIRE(vptr->isSimpleVariantLocus() );
    BOOST_REQUIRE( vptr->is_indel() );
    BOOST_REQUIRE( ! vptr->is_snv() );
    BOOST_REQUIRE( vptr->is_normalized());
    BOOST_REQUIRE( vptr->is_match_reference());
    BOOST_REQUIRE_EQUAL(vptr->pos, 21);
    BOOST_REQUIRE_EQUAL(vptr->alt.size(),1u);
    BOOST_REQUIRE_EQUAL(vptr->alt[0],"AGCT");

    BOOST_REQUIRE( vcfs.next() );
    vptr = vcfs.get_record_ptr();
    assert(vptr != nullptr);

    // VCF record test for complex val
    // chrQ   29    INVALID_ALT    A    <DEL>    ...
    // this should not be parsed as an SNV or indel
    BOOST_REQUIRE(not vptr->isSimpleVariantLocus());
    BOOST_REQUIRE_THROW( vcfs.next(), illumina::common::LogicException );

    // VCF record test for
    // chrQ    30    MISMATCH_REF_SNP    T    A    ...
    // testing that mismatched ref is caught
    BOOST_REQUIRE( vptr->isSimpleVariantLocus() );
    BOOST_REQUIRE( ! vptr->is_indel() );
    BOOST_REQUIRE( vptr->is_snv() );
    BOOST_REQUIRE( vptr->is_normalized());
    BOOST_REQUIRE( ! vptr->is_match_reference());
    BOOST_REQUIRE_EQUAL(vptr->pos, 30);
    BOOST_REQUIRE_EQUAL(vptr->alt.size(),1u);
    BOOST_REQUIRE_EQUAL(vptr->alt[0],"A");

    BOOST_REQUIRE_THROW( vcfs.next(), illumina::common::LogicException );
    vptr = vcfs.get_record_ptr();
    assert(vptr != nullptr);

    // VCF record test for
    // chrQ    32    MISMATCH_REF_INDEL  GTA    G    ...
    // testing that mismatched ref is caught
    BOOST_REQUIRE(vptr->isSimpleVariantLocus() );
    BOOST_REQUIRE( vptr->is_indel() );
    BOOST_REQUIRE( ! vptr->is_snv() );
    BOOST_REQUIRE( vptr->is_normalized());
    BOOST_REQUIRE( ! vptr->is_match_reference());
    BOOST_REQUIRE_EQUAL(vptr->pos, 32);
    BOOST_REQUIRE_EQUAL(vptr->alt.size(),1u);
    BOOST_REQUIRE_EQUAL(vptr->alt[0],"G");

    // testing that next returns false after last record
    BOOST_REQUIRE( ! vcfs.next() );
    vptr = vcfs.get_record_ptr();
    assert(vptr == nullptr);
}


BOOST_AUTO_TEST_SUITE_END()

