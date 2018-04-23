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
#include "boost/algorithm/string.hpp"

#include "VariantPhaser.hh"

struct DummyVariantSink : public variant_pipe_stage_base
{
    DummyVariantSink() : variant_pipe_stage_base() {}

    void process(std::unique_ptr<GermlineSiteLocusInfo> siteLocus) override
    {
        _siteLoci.push_back(std::move(siteLocus));
    }
    void process(std::unique_ptr<GermlineIndelLocusInfo> indelLocus) override
    {
        _indelLoci.push_back(std::move(indelLocus));
    }

    bool check(const unsigned sampleIndex, const bool isSnv, const pos_t pos,
               const bool isPhased, const uint8_t allele0Index, const uint8_t allele1Index)
    {
        if (isSnv)
        {
            for (const auto& locus : _siteLoci)
            {
                if (locus->pos != pos) continue;

                const auto& maxGenotype(locus->getSample(sampleIndex).max_gt());
                return (maxGenotype.getIsPhased() == isPhased)
                       && (maxGenotype.getAllele0Index() == allele0Index)
                       && (maxGenotype.getAllele1Index() == allele1Index);
            }
        }
        else
        {
            for (const auto& locus : _indelLoci)
            {
                if (locus->pos != pos) continue;

                const auto& maxGenotype(locus->getSample(sampleIndex).max_gt());
                return (maxGenotype.getIsPhased() == isPhased)
                       && (maxGenotype.getAllele0Index() == allele0Index)
                       && (maxGenotype.getAllele1Index() == allele1Index);
            }
        }

        return false;
    }

    std::vector<std::unique_ptr<GermlineSiteLocusInfo>> _siteLoci;
    std::vector<std::unique_ptr<GermlineIndelLocusInfo>> _indelLoci;
};



static
std::unique_ptr<GermlineDiploidIndelLocusInfo>
getDeletionLocusInfo(
    const starling_deriv_options& dopt,
    const unsigned sampleCount,
    const unsigned sampleIndex,
    const ActiveRegionId activeRegionId,
    const pos_t pos,
    const unsigned indelLength,
    const uint8_t haplotypeId,
    bool isHom
)
{
    IndelKey indelKey;
    indelKey.pos=pos;
    indelKey.type=INDEL::INDEL;
    indelKey.deletionLength=indelLength;

    IndelData indelData(sampleCount,indelKey);
    indelData.getSampleData(sampleIndex).haplotypeId = haplotypeId;

    std::unique_ptr<GermlineDiploidIndelLocusInfo> indelInfo(new GermlineDiploidIndelLocusInfo(dopt.gvcf, sampleCount));
    indelInfo->addAltIndelAllele(indelKey, indelData);

    auto& indelSampleInfo(indelInfo->getSample(sampleIndex));
    indelSampleInfo.setActiveRegionId(activeRegionId);
    if (!isHom)
        indelSampleInfo.max_gt().setGenotypeFromAlleleIndices(0,1);
    else
        indelSampleInfo.max_gt().setGenotypeFromAlleleIndices(1,1);

    indelInfo->getSample(sampleIndex).max_gt().setAllele1HaplotypeId(haplotypeId);
    return indelInfo;
}



static
std::unique_ptr<GermlineDiploidSiteLocusInfo>
getSnvLocusInfo(
    const starling_deriv_options& dopt,
    const unsigned sampleCount,
    const unsigned sampleIndex,
    const ActiveRegionId activeRegionId,
    const pos_t pos,
    const char refBaseChar,
    const char altBaseChar,
    const uint8_t allele0HaplotypeId,
    const uint8_t allele1HaplotypeId,
    bool isHom)
{
    std::unique_ptr<GermlineDiploidSiteLocusInfo> siteInfo(
        new GermlineDiploidSiteLocusInfo(dopt.gvcf, sampleCount, pos, base_to_id(refBaseChar)));

    siteInfo->addAltSiteAllele(static_cast<BASE_ID::index_t>(base_to_id(altBaseChar)));

    auto& siteSampleInfo(siteInfo->getSample(sampleIndex));
    siteSampleInfo.setActiveRegionId(activeRegionId);
    auto& maxGenotype(siteSampleInfo.max_gt());
    if (isHom)
        maxGenotype.setGenotypeFromAlleleIndices(1, 1);
    else
        maxGenotype.setGenotypeFromAlleleIndices(0, 1);

    maxGenotype.setAllele0HaplotypeId(allele0HaplotypeId);
    maxGenotype.setAllele1HaplotypeId(allele1HaplotypeId);

    return siteInfo;
}


static
starling_options
getMockOptions()
{
    starling_options opt;
    opt.alignFileOpt.alignmentFilenames.push_back("sample.bam");
    opt.isUseVariantPhaser = true;
    return opt;
}


BOOST_AUTO_TEST_SUITE( variantPhaserTestSuite )

BOOST_AUTO_TEST_CASE( multiSamplePhasingTest )
{
    const starling_options opt = getMockOptions();
    const starling_deriv_options dopt(opt);

    std::shared_ptr<DummyVariantSink> next(new DummyVariantSink);

    const unsigned sampleCount(2);
    const ActiveRegionId activeRegionId(1);

    VariantPhaser phaser(opt, sampleCount, next);

    // =============== Sample 1 ==================
    // hap0 (ref):  CAAACAAAAAAACAAAAAAAAACAAAAAAT
    // hap1:        CAAA--AAAAAAGCAAAAAAAACAAAAAAC
    // hap2:        CAAACAAAAAAAC-AAAAAAAACAAAAAAC
    // =============== Sample 2 ==================
    // hap0 (ref):  CAAACAAAAAAACAAAAAAAAACAAAAAAT
    // hap1:        C--ACAAAAAAAGCAAAAAAAACAAAAAAT
    // hap2:        CAAACAAAAATACCAAAAAAAACAAAAAAT

    const unsigned sampleIndex0(0);
    const unsigned sampleIndex1(1);

    // sample: 2, pos: 1, deletionLength: 2, hap 1, genotype: 0|1
    auto indelInfo1(getDeletionLocusInfo(dopt, sampleCount, sampleIndex1, activeRegionId, 1, 0, 1, false));
    phaser.process(std::move(indelInfo1));

    // sample: 1, pos: 4, deletionLength: 2, hap 1, genotype: 0|1
    auto indelInfo2(getDeletionLocusInfo(dopt, sampleCount, sampleIndex0, activeRegionId, 4, 0, 1, false));
    phaser.process(std::move(indelInfo2));

    // sample: 2, pos: 10, A->T, hap 2, genotype: 1|0
    auto snvInfo2(getSnvLocusInfo(dopt, sampleCount, sampleIndex1, activeRegionId, 10, 'A', 'T', 0, 2, false));
    phaser.process(std::move(snvInfo2));

    // sample: 1, pos: 12, C->G, hap 1, genotype: 0|1
    auto snvInfo1(getSnvLocusInfo(dopt, sampleCount, sampleIndex0, activeRegionId, 12, 'C', 'G', 0, 1, false));
    phaser.process(std::move(snvInfo1));

    // sample: 1, pos: 13, deletionLength: 1, hap 2, genotype: 1|0
    auto indelInfo3(getDeletionLocusInfo(dopt, sampleCount, sampleIndex0, activeRegionId, 13, 0, 2, false));
    phaser.process(std::move(indelInfo3));

    // sample: 1, pos: 13, A->C, hap 1, genotype: 0|1
    auto snvInfo3(getSnvLocusInfo(dopt, sampleCount, sampleIndex0, activeRegionId, 13, 'A', 'C', 0, 1, false));
    phaser.process(std::move(snvInfo3));

    // sample: 1, pos: 29, T->C, hap 1 and 2, genotype: 1/1
    // alt1HaplotypeId==3 means that T->C appears both hap1 and hap2
    auto snvInfo4(getSnvLocusInfo(dopt, sampleCount, sampleIndex0, activeRegionId, 29, 'T', 'C', 0, 3, true));
    phaser.process(std::move(snvInfo4));

    phaser.flush();

    // sample: 2, indel pos: 1 deletionLength: 2, hap 1, genotype: 0|1
    BOOST_REQUIRE(next->check(sampleIndex1, false, 1, true, 0, 1));
    // sample: 1, indel pos: 4, deletionLength: 2, hap 1, genotype: 0|1
    BOOST_REQUIRE(next->check(sampleIndex0, false, 4, true, 0, 1));
    // sample: 2, pos: 10, A->T, hap 2, genotype: 1|0
    BOOST_REQUIRE(next->check(sampleIndex1, true, 10, true, 1, 0));
    // sample: 1, pos: 12, C->G, hap 1, genotype: 0|1
    BOOST_REQUIRE(next->check(sampleIndex0, true, 12, true, 0, 1));
    // sample: 1, pos: 13, deletionLength: 1, hap 2, genotype: 1|0
    BOOST_REQUIRE(next->check(sampleIndex0, false, 13, true, 1, 0));
    // sample: 1, pos: 13, A->C, hap 1, genotype: 0|1
    BOOST_REQUIRE(next->check(sampleIndex0, true, 13, true, 0, 1));
    // sample: 1, pos: 29, T->C, hap 1 and 2, genotype: 1/1
    BOOST_REQUIRE(next->check(sampleIndex0, true, 29, false, 1, 1));
}

BOOST_AUTO_TEST_CASE( phasingConflictTest )
{
    const starling_options opt = getMockOptions();
    const starling_deriv_options dopt(opt);

//        std::shared_ptr<variant_pipe_stage_base> next(new DummyVariantSink);
    std::shared_ptr<DummyVariantSink> next(new DummyVariantSink);

    const unsigned sampleCount(1);
    const unsigned sampleIndex(0);
    const ActiveRegionId activeRegionId(1);

    VariantPhaser phaser(opt, sampleCount, next);

    // hap0 (ref):  CAAATCAT
    // hap1:        CTAA-CAT
    // hap1:        CAAA-CGT

    // pos: 1, A->T, hap 1, genotype: 0|1
    auto snvInfo1(getSnvLocusInfo(dopt, sampleCount, sampleIndex, activeRegionId, 1, 'A', 'T', 0, 1, false));
    phaser.process(std::move(snvInfo1));

    // pos: 2, A->T, not in selected haplotype, genotype should be 0/1 (unphased)
    auto conflictSnvInfo(getSnvLocusInfo(dopt, sampleCount, sampleIndex, activeRegionId, 2, 'A', 'T', 0, 0, false));
    phaser.process(std::move(conflictSnvInfo));

    // pos: 4, deletionLength: 1, hap 1 and 2 (id==3), genotype should be 0/1 (unphased)
    // hom variant in active region, but called het
    auto conflictIndelInfo(getDeletionLocusInfo(dopt, sampleCount, sampleIndex, activeRegionId, 4, 0, 3, false));
    phaser.process(std::move(conflictIndelInfo));

    // pos: 6, A->G, hap 2, genotype: 1|0
    auto snvInfo2(getSnvLocusInfo(dopt, sampleCount, sampleIndex, activeRegionId, 6, 'A', 'G', 0, 2, false));
    phaser.process(std::move(snvInfo2));

    phaser.flush();

    // pos: 1, A->T, hap 1, genotype: 0|1
    BOOST_REQUIRE(next->check(sampleIndex, true, 1, true, 0, 1));
    // pos: 2, A->T, not in selected haplotype, genotype should be 0/1 (unphased)
    BOOST_REQUIRE(next->check(sampleIndex, true, 2, false, 0, 1));
    // pos: 4, deletionLength: 1, hap 1 and 2 (id==3), genotype should be 0/1 (unphased)
    BOOST_REQUIRE(next->check(sampleIndex, false, 4, false, 0, 1));
    // pos: 6, A->G, hap 2, genotype: 1|0
    BOOST_REQUIRE(next->check(sampleIndex, true, 6, true, 1, 0));
}

BOOST_AUTO_TEST_SUITE_END()
