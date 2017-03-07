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

#include <applications/starling/LocusReportInfoUtil.hh>
#include "boost/test/unit_test.hpp"
#include "boost/algorithm/string.hpp"

#include "test_config.h"

#include "gvcf_writer.hh"
#include "ScoringModelManager.hh"

#include "starling_common/pos_basecall_buffer.hh"

BOOST_AUTO_TEST_SUITE( gvcf_writer_test )

BOOST_AUTO_TEST_CASE( test_removeCommonPrefix )
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

#if 0
static void
insert_read(
    const char* read,
    pos_t position,
    pos_basecall_buffer& bc_buff)
{
    pos_t insert_pos = position;
    while (*read)
    {
        base_call bc(base_to_id(*read), 30, true, 0, 0, false, false, false,
                     insert_pos == position, !*(read + 1));
        bc_buff.insert_pos_basecall(insert_pos, true, bc);
        insert_pos++;
        read++;
    }

}




BOOST_AUTO_TEST_CASE( unphased_flag_written )
{
    reference_contig_segment rcs;
    rcs.seq() = "ACGTACGTACGT";
    pos_basecall_buffer bc_buff(rcs);

    auto r1 = "ACGTACGTAC";
    auto r2 = "AGGTACGTAC";
    pos_t read_pos = 0;
    pos_t snp_pos = 1;

    // add 2 haplotypes of reads
    for (int i = 0; i < 10; i++)
        insert_read(r1, read_pos, bc_buff);
    for (int i = 0; i < 10; i++)
        insert_read(r2, read_pos, bc_buff);

    starling_options opt;
    opt.gvcf.is_skip_header = true;
    opt.is_user_genome_size = true;
    opt.user_genome_size = rcs.seq().size();
    opt.bam_seq_name = "dummy";
    starling_deriv_options dopt(opt,rcs);

    const unsigned sampleCount(1);
    const unsigned sampleIndex(0);
    //const snp_pos_info& spi(bc_buff.get_pos(snp_pos));
    std::unique_ptr<GermlineDiploidSiteLocusInfo> locusPtr(new GermlineDiploidSiteLocusInfo(dopt.gvcf, sampleCount, snp_pos, base_to_id(rcs.get_base(snp_pos))));

    auto& siteSampleInfo(locusPtr->getSiteSample(sampleIndex));
    siteSampleInfo.max_gt = DIGT::get_gt_with_alleles(base_to_id(r1[snp_pos]),base_to_id(r2[snp_pos]));
    locusPtr->anyVariantAlleleQuality = 60;
    locusPtr->addAltSiteAllele(BASE_ID::G);

    RegionTracker regions;
    const std::vector<std::string> sampleNames = { "SAMPLE" };
    std::stringstream os;

    gvcf_deriv_options gvcf_options(dopt.gvcf);
    gvcf_options.chrom_depth["dummy"] = 30.734;
    ScoringModelManager cm(opt, gvcf_options);


    gvcf_writer writer(opt, dopt, rcs, regions, sampleNames, &os, cm);
    writer.process(std::move(locusPtr));

    std::string x = os.str();
    BOOST_REQUIRE(x.length() > 0);
    std::vector<std::string> lines;
    boost::split(lines, x, boost::is_any_of("\n"));
    BOOST_REQUIRE(lines.size() == 3);

    std::vector<std::string> strs;
    boost::split(strs, lines[1], boost::is_any_of("\t"));
    BOOST_REQUIRE(strs.size() > 7);
}

#endif
