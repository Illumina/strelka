// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
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

#include "codon_phaser.hh"



static void insert_read(const char* read, pos_t position,
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



BOOST_AUTO_TEST_SUITE( codon_phaser )

class dummy_variant_sink : public variant_pipe_stage_base
{
public:
    dummy_variant_sink() : variant_pipe_stage_base() {}
    std::vector<std::unique_ptr<GermlineDiploidSiteLocusInfo>> the_sites;
    std::vector<std::unique_ptr<GermlineIndelLocusInfo>> the_indels;
    void process(std::unique_ptr<GermlineSiteLocusInfo> site) override
    {
        auto si(downcast<GermlineDiploidSiteLocusInfo>(std::move(site)));
        if (si->is_het() || si->is_hetalt() ) the_sites.push_back(std::move(si));
    }
    void process(std::unique_ptr<GermlineIndelLocusInfo> ii) override
    {
        the_indels.push_back(std::move(ii));
    }
};



static
std::shared_ptr<dummy_variant_sink>
getTestSink(
    const char* refSeq,
    const char* read1Seq,
    const char* read2Seq,
    const int phasingWindow,
    const int depth)
{
    assert(strlen(read1Seq) == strlen(read2Seq));

    reference_contig_segment rcs;
    rcs.seq() = refSeq;
    pos_basecall_buffer bc_buff(rcs);

    pos_t read_pos = 0;

    // add 2 haplotypes of reads
    for (int i = 0; i < depth; i++)
    {
        insert_read(read1Seq, read_pos, bc_buff);
        insert_read(read2Seq, read_pos, bc_buff);
    }
    starling_options opt;
    opt.phasing_window = phasingWindow;
    opt.do_codon_phasing = true;
    opt.is_user_genome_size = true;
    opt.user_genome_size = rcs.seq().size();
    opt.bam_seq_name = "dummy";
    starling_deriv_options dopt(opt, rcs);

    std::shared_ptr<dummy_variant_sink> next(new dummy_variant_sink);
    Codon_phaser phaser(opt, bc_buff, rcs, std::dynamic_pointer_cast<variant_pipe_stage_base>(next));

    const unsigned sampleCount(1);

    for (int i = 0; read1Seq[i]; i++)
    {
        const snp_pos_info& spi(bc_buff.get_pos(read_pos + i));
        std::unique_ptr<GermlineDiploidSiteLocusInfo> si(
            new GermlineDiploidSiteLocusInfo(dopt.gvcf, sampleCount, read_pos + i, rcs.get_base(read_pos + i), spi, 30));
        si->allele.is_covered = si->allele.is_used_covered = true;
        si->dgt.ref_gt = base_to_id(si->ref);

        si->allele.max_gt = DIGT::get_gt_with_alleles(base_to_id(read1Seq[i]), base_to_id(read2Seq[i]));
        si->dgt.is_snp = si->ref != read1Seq[i] || si->ref != read2Seq[i];

        phaser.process(std::move(si));
    }
    phaser.flush();

    return next;
}

static
void
positiveTest(
    const char* refSeq,
    const char* read1Seq,
    const char* read2Seq,
    const int phasingWindow,
    const char* expectedPhasedAlt)
{
    const auto next = getTestSink(refSeq, read1Seq, read2Seq, phasingWindow, 10);
    BOOST_CHECK_EQUAL(expectedPhasedAlt, next->the_sites.front()->phased_alt);
}

// positive tests
BOOST_AUTO_TEST_CASE( simple_3mer )
{
    positiveTest("ACGTACGTACGT", "ACGTACGTAC", "ACGTGCTTAC", 3, "GCT");
}

BOOST_AUTO_TEST_CASE( two_adjacent_3mers )
{
    positiveTest("ACGTACGTACGTACGT", "ACGTACGTACGTACGT", "ACGTGCTTGCTTACGT", 3, "GCTTGCT");
}

BOOST_AUTO_TEST_CASE( handles_snps_at_start )
{
    positiveTest("ACGTACGTACGT", "ACGTACGTAC", "GGGTACGTAC", 3, "GG");
}


BOOST_AUTO_TEST_CASE( respects_phasing_window )
{
    positiveTest("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGTACGT",
                 "ACGTAAGAACGGAGGTACGTACGT", 5, "AGAACGGAG");
}

BOOST_AUTO_TEST_CASE(test_overlapping_phased_snps)
{
    positiveTest("ACGTACGTACGT", "ACATGCGTAC", "ACCTCCGTAC", 3, "ATG,CTC");
}

BOOST_AUTO_TEST_CASE(test_overlapping_phased_snps_different_size)
{
    positiveTest("ACGTACGTACGT", "ACGAAAGTAC", "ACATCCGTAC", 5, "ATCC,GAAA");
}


// negative tests
BOOST_AUTO_TEST_CASE( just_one_snp )
{
    const auto next = getTestSink("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", "ACGTACGT", "ACGGACGT", 3, 10);

    for (auto& phased_variant : next->the_sites)
    {
        BOOST_CHECK(!phased_variant->filters.any());
        BOOST_CHECK(!phased_variant->allele.is_phased_region);
    }
    BOOST_CHECK_EQUAL(next->the_sites.size(), 1);
}

BOOST_AUTO_TEST_CASE( read_break_causes_phasing_conflict )
{
    reference_contig_segment rcs;
    rcs.seq() = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    pos_basecall_buffer bc_buff(rcs);

    auto r1 = "ACGTACGTACGTACGT";
    auto r1_nosnp1 = "ACGTACGT";
    auto r1_nosnp2 = "ACGTACGT";
    auto r2 = "ACGTACGGAGGTACGT";

    pos_t read_pos = 0;
    pos_t read_pos2 = strlen(r1) + read_pos;
    pos_t read_pos_snp = read_pos;

    // add 2 haplotypes of reads
    for (int i = 0; i < 10; i++)
    {
        insert_read(r1_nosnp1, read_pos, bc_buff);
        insert_read(r1_nosnp2, read_pos2, bc_buff);
    }

    for (int i = 0; i < 10; i++)
    {
        insert_read(r2, read_pos_snp, bc_buff);
    }

    starling_options opt;
    opt.phasing_window = 3;
    opt.do_codon_phasing = true;
    opt.is_user_genome_size = true;
    opt.user_genome_size = rcs.seq().size();
    opt.bam_seq_name = "dummy";
    starling_deriv_options dopt(opt, rcs);

    std::shared_ptr<dummy_variant_sink> next(new dummy_variant_sink);
    Codon_phaser phaser(opt, bc_buff, rcs, std::dynamic_pointer_cast<variant_pipe_stage_base>(next));
    dummy_variant_sink& sink(*next);

    /// TODO STREL-125 generalize to mulit-sample or replace codon phaser
    const unsigned sampleCount(1);
    const unsigned sampleIndex(0);

    for (int i = 0; r1[i]; i++)
    {
        const snp_pos_info& spi(bc_buff.get_pos(read_pos + i));
        std::unique_ptr<GermlineDiploidSiteLocusInfo> si(new GermlineDiploidSiteLocusInfo(dopt.gvcf, sampleCount, read_pos + i, rcs.get_base(read_pos + i), spi, 30));
        auto& sampleInfo(si->getSample(sampleIndex));
        si->allele.is_covered = si->allele.is_used_covered = true;
        sampleInfo.genotypeQualityPolymorphic = si->dgt.genome.snp_qphred = sampleInfo.empiricalVariantScore = 40;
        si->dgt.ref_gt = base_to_id(si->ref);

        si->allele.max_gt = DIGT::get_gt_with_alleles(base_to_id(r1[i]),base_to_id(r2[i]));
        si->dgt.is_snp = si->ref != r1[i] || si->ref != r2[i];

        phaser.process(std::move(si));
    }
    phaser.flush();
    for (auto& site : sink.the_sites)
    {
        BOOST_CHECK(! site->is_het() || site->filters.test(GERMLINE_VARIANT_VCF_FILTERS::PhasingConflict));
        BOOST_CHECK(!site->allele.is_phased_region);
    }
}

BOOST_AUTO_TEST_CASE( low_depth_doesnt_phase )
{
    const auto next = getTestSink("ACGTACGTACGT", "ACGTACGTAC", "ACGTGCTTAC", 3, 4);

    BOOST_CHECK_EQUAL(2, next->the_sites.size());
    BOOST_CHECK(next->the_sites.front()->allele.is_phasing_insufficient_depth);
    BOOST_CHECK(next->the_sites.back()->allele.is_phasing_insufficient_depth);

}



// TODO: write tests for:
// allele imbalanced phasing
// nonsense phasing (i.e. the called SNPs are not in the top most common alleles)

BOOST_AUTO_TEST_SUITE_END()
