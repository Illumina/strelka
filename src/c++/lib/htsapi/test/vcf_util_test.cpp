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

#include "htsapi/vcf_util.hh"

#include "boost/test/unit_test.hpp"


static
void
testGtToAllele(
    const unsigned genotypeIndex,
    const unsigned expectedAllele0Index,
    const unsigned expectedAllele1Index)
{
    uint8_t a0,a1;
    VcfGenotypeUtil::getAlleleIndices(genotypeIndex,a0,a1);
    BOOST_REQUIRE_EQUAL(static_cast<unsigned>(a0), expectedAllele0Index);
    BOOST_REQUIRE_EQUAL(static_cast<unsigned>(a1), expectedAllele1Index);
}


static
void
testAlleleToGt(
    const unsigned expectedGenotypeIndex,
    const unsigned allele0Index,
    const unsigned allele1Index)
{
    const unsigned genotypeIndex = VcfGenotypeUtil::getGenotypeIndex(allele0Index, allele1Index);
    BOOST_REQUIRE_EQUAL(genotypeIndex, expectedGenotypeIndex);
}



BOOST_AUTO_TEST_SUITE( test_vcf_util )


BOOST_AUTO_TEST_CASE( vcf_genotype_to_allele_test )
{
    testGtToAllele(0, 0, 0);
    testGtToAllele(1, 0, 1);
    testGtToAllele(2, 1, 1);
    testGtToAllele(3, 0, 2);
    testGtToAllele(4, 1, 2);
    testGtToAllele(5, 2, 2);
    testGtToAllele(6, 0, 3);
}

BOOST_AUTO_TEST_CASE( vcf_allele_to_genotype_test )
{
    testAlleleToGt(0, 0, 0);
    testAlleleToGt(1, 0, 1);
    testAlleleToGt(2, 1, 1);
    testAlleleToGt(3, 0, 2);
    testAlleleToGt(4, 1, 2);
    testAlleleToGt(5, 2, 2);
    testAlleleToGt(6, 0, 3);
}


BOOST_AUTO_TEST_SUITE_END()

