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

#include "testConfig.h"

#include "starling_pos_processor_util.hh"


static
std::string
getTestpath()
{
    static const std::string testPath(std::string(TEST_DATA_PATH) + "/subregionTest.bed.gz");
    return testPath;
}


BOOST_AUTO_TEST_SUITE( sppr_util_test_suite )

BOOST_AUTO_TEST_CASE( vcf_to_allele_test )
{
    static const unsigned max_indel_size(50);
    IndelObservation obs;
    bool isConverted;

    {
        // "chr20 4329513 . ATT AT,A 9114 PASS LEN=1,2;TYPE=del,del GT 2|1"
        vcf_record vr;
        vr.chrom = "chr20";
        vr.pos = 4329513;
        vr.ref = "ATT";
        vr.alt.push_back("AT");
        vr.alt.push_back("A");

        isConverted = convert_vcfrecord_to_indel_allele(max_indel_size,vr,0,obs);
        IndelKey k0expect(4329513,INDEL::INDEL,1);
        BOOST_REQUIRE(isConverted);
        BOOST_REQUIRE_EQUAL(obs.key,k0expect);

        isConverted = convert_vcfrecord_to_indel_allele(max_indel_size,vr,1,obs);
        IndelKey k1expect(4329513,INDEL::INDEL,2);
        BOOST_REQUIRE(isConverted);
        BOOST_REQUIRE_EQUAL(obs.key,k1expect);
    }

    {
        // "chr20 4329513 . ATT A,ATTT,AG"
        vcf_record vr;
        vr.chrom = "chr20";
        vr.pos = 4329513;
        vr.ref = "ATT";
        vr.alt.push_back("A");
        vr.alt.push_back("ATTT");
        vr.alt.push_back("AG");

        isConverted = convert_vcfrecord_to_indel_allele(max_indel_size,vr,0,obs);
        IndelKey k0expect(4329513,INDEL::INDEL,2);
        BOOST_REQUIRE(isConverted);
        BOOST_REQUIRE_EQUAL(obs.key,k0expect);

        isConverted = convert_vcfrecord_to_indel_allele(max_indel_size,vr,1,obs);
        IndelKey k1expect(4329513,INDEL::INDEL, 0, "T");
        BOOST_REQUIRE(isConverted);
        BOOST_REQUIRE_EQUAL(obs.key,k1expect);

        isConverted = convert_vcfrecord_to_indel_allele(max_indel_size,vr,2,obs);
        IndelKey k2expect(4329513,INDEL::INDEL, 2, "G");
        BOOST_REQUIRE(isConverted);
        BOOST_REQUIRE_EQUAL(obs.key,k2expect);
    }
}


BOOST_AUTO_TEST_CASE( subregions_from_bed_test )
{
    static const std::string regionChrom("chrTest");
    static const known_pos_range2 regionRange(10000,30000);
    std::vector<known_pos_range2> subRegionRanges;
    getSubRegionsFromBedTrack(getTestpath(), regionChrom, regionRange, subRegionRanges);

    BOOST_REQUIRE_EQUAL(subRegionRanges.size(),3);
    BOOST_REQUIRE_EQUAL(subRegionRanges[0],known_pos_range2(10000,10200));
    BOOST_REQUIRE_EQUAL(subRegionRanges[1],known_pos_range2(18000,22000));
    BOOST_REQUIRE_EQUAL(subRegionRanges[2],known_pos_range2(29800,30000));
}


BOOST_AUTO_TEST_SUITE_END()
