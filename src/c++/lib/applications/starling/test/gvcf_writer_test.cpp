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

#include "gvcf_writer.hh"
#include "LocusReportInfoUtil.hh"


BOOST_AUTO_TEST_SUITE( gvcf_writer_test_suite )

BOOST_AUTO_TEST_CASE( removeCommonPrefixTest )
{
    reference_contig_segment ref;
    ref.seq() = "TAAGTGAAGTATTTTTTTTTTTTT";

    auto indelKey1 = IndelKey(5, INDEL::INDEL, 11, "A");
    auto indelKey2 = IndelKey(4, INDEL::INDEL, 9);
    auto indelKey3 = IndelKey(1, INDEL::INDEL, 5);

    const unsigned sampleCount(1);

    std::vector<GermlineIndelAlleleInfo> indelAlleles;
    indelAlleles.push_back(GermlineIndelAlleleInfo(indelKey1, IndelData(sampleCount, indelKey1)));
    indelAlleles.push_back(GermlineIndelAlleleInfo(indelKey2, IndelData(sampleCount, indelKey1)));
    indelAlleles.push_back(GermlineIndelAlleleInfo(indelKey3, IndelData(sampleCount, indelKey1)));

    const unsigned commonPrefixLength(4);

    OrthogonalAlleleSetLocusReportInfo locusReportInfo;

    getLocusReportInfoFromAlleles(ref, indelAlleles, commonPrefixLength, locusReportInfo);

    BOOST_REQUIRE_EQUAL(locusReportInfo.vcfPos, 1+commonPrefixLength);
    BOOST_REQUIRE_EQUAL(ALIGNPATH::apath_to_cigar(locusReportInfo.altAlleles[0].vcfCigar), "1M11D1I");
    BOOST_REQUIRE_EQUAL(ALIGNPATH::apath_to_cigar(locusReportInfo.altAlleles[1].vcfCigar), "1M9D2M");
    BOOST_REQUIRE_EQUAL(ALIGNPATH::apath_to_cigar(locusReportInfo.altAlleles[2].vcfCigar), "1M5D6M");
}

BOOST_AUTO_TEST_SUITE_END()
