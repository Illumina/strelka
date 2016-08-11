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

#include "indel_overlapper.hh"

#include "ScoringModelManager.hh"


struct dummy_variant_sink : public variant_pipe_stage_base
{
    dummy_variant_sink() : variant_pipe_stage_base() {}

    void process(std::unique_ptr<GermlineSiteLocusInfo> si) override
    {
        the_sites.push_back(std::move(si));
    }
    void process(std::unique_ptr<GermlineIndelLocusInfo> ii) override
    {
        the_indels.push_back(std::move(ii));
    }

    std::vector<std::unique_ptr<GermlineSiteLocusInfo>> the_sites;
    std::vector<std::unique_ptr<GermlineIndelLocusInfo>> the_indels;
};



BOOST_AUTO_TEST_SUITE( indel_overlapper_test )

/// add a single indel to the indel_overlapper, this is
/// a simple demo to unit test the overlapper
///
BOOST_AUTO_TEST_CASE( simple_indel_test )
{
    // fake various high-level data structures with as many defaults as possible
    reference_contig_segment rcs;
    rcs.seq() = "ACGGGGTTGGACGATGCTACGATCGATCGCGTACCTACGATCGACTACGACTGCGACGATCGACGATCGACGATCGATCGATCGACGTACGACACGTACGATCGATCGATCGATCGACTCGATCAGCTCATGCATCG";

    starling_options opt;
    opt.bam_seq_name="chr1";
    opt.is_user_genome_size = true;
    opt.user_genome_size = rcs.seq().size();

    starling_deriv_options dopt(opt,rcs);

    ScoringModelManager cm(opt, dopt.gvcf);

    std::shared_ptr<variant_pipe_stage_base> next(new dummy_variant_sink);
    indel_overlapper overlap(cm, rcs, next);

    IndelKey indelKey;
    const IndelData indelData(1,indelKey);
    const GermlineDiploidIndelSimpleGenotypeInfoCore dindel;
    const AlleleReportInfo indelReportInfo;

    indelKey.pos=6;
    indelKey.type=INDEL::INDEL;
    indelKey.deletionLength=2;

    const unsigned sampleCount(1);
    std::unique_ptr<GermlineDiploidIndelLocusInfo> ii(new GermlineDiploidIndelLocusInfo(dopt.gvcf, sampleCount, indelKey,indelData,dindel,indelReportInfo));
    overlap.process(std::move(ii));

    overlap.flush();
    // just testing that we don't assert for now...
}


/// Reduced recreation of STARKA-393 failure:
BOOST_AUTO_TEST_CASE( conflicting_indel_test )
{
    // fake various high-level data structures with as many defaults as possible
    reference_contig_segment rcs;
    rcs.seq() = "ACGGGGTTGGACGATGCTACGATCGATCGCGTACCTACGATCGACTACGACTGCGACGATCGACGATCGACGATCGATCGATCGACGTACGACACGTACGATCGATCGATCGATCGACTCGATCAGCTCATGCATCG";

    starling_options opt;
    opt.bam_seq_name="chr1";
    opt.is_user_genome_size = true;
    opt.user_genome_size = rcs.seq().size();

    starling_deriv_options dopt(opt,rcs);

    ScoringModelManager cm(opt, dopt.gvcf);

    std::shared_ptr<variant_pipe_stage_base> next(new dummy_variant_sink);
    indel_overlapper overlap(cm, rcs, next);

    IndelKey iks[] =
    {
        IndelKey(10,INDEL::INDEL,10),
        IndelKey(15,INDEL::INDEL,30),
        IndelKey(20,INDEL::INDEL,1),
        IndelKey(25,INDEL::INDEL,1),
    };

    int max_gts[] = { 2,0,2,2 };

    IndelKey ik;
    IndelData indelData(1,ik);
    indelData.isForcedOutput = true;
    GermlineDiploidIndelSimpleGenotypeInfoCore dindel;
    const AlleleReportInfo indelReportInfo;

    const unsigned sampleCount(1);
    const unsigned keyCount(sizeof(iks)/sizeof(IndelKey));
    for (unsigned i(0); i<keyCount; ++i)
    {
        ik=iks[i];
        dindel.max_gt=max_gts[i];
        dindel.max_gt_poly=max_gts[i];
        std::unique_ptr<GermlineDiploidIndelLocusInfo> ii(new GermlineDiploidIndelLocusInfo(dopt.gvcf, sampleCount, ik,indelData,dindel,indelReportInfo));
        overlap.process(std::move(ii));
    }

    overlap.flush();
    // just testing that we don't assert for now...
}


/// Additional test added in light of STARKA-393 failure:
BOOST_AUTO_TEST_CASE( conflicting_indel_test2 )
{
    // fake various high-level data structures with as many defaults as possible
    reference_contig_segment rcs;
    rcs.seq() = "ACGGGGTTGGACGATGCTACGATCGATCGCGTACCTACGATCGACTACGACTGCGACGATCGACGATCGACGATCGATCGATCGACGTACGACACGTACGATCGATCGATCGATCGACTCGATCAGCTCATGCATCG";

    starling_options opt;
    opt.bam_seq_name="chr1";
    opt.is_user_genome_size = true;
    opt.user_genome_size = rcs.seq().size();

    starling_deriv_options dopt(opt,rcs);

    ScoringModelManager cm(opt, dopt.gvcf);

    std::shared_ptr<dummy_variant_sink> next(new dummy_variant_sink);
    indel_overlapper overlap(cm, rcs, next);

    IndelKey indelKeys[] =
    {
        IndelKey(10,INDEL::INDEL,10),
        IndelKey(12,INDEL::INDEL,1),
        IndelKey(15,INDEL::INDEL,1),
        IndelKey(18,INDEL::INDEL,1),
    };

    int max_gts[] = { 2,2,0,2 };

    IndelKey indelKey;
    IndelData indelData(1,indelKey);
    indelData.isForcedOutput = true;
    GermlineDiploidIndelSimpleGenotypeInfoCore dindel;
    const AlleleReportInfo indelReportInfo;

    const unsigned sampleCount(1);
    const unsigned keyCount(sizeof(indelKeys)/sizeof(IndelKey));
    for (unsigned i(0); i<keyCount; ++i)
    {
        indelKey=indelKeys[i];
        dindel.max_gt=max_gts[i];
        dindel.max_gt_poly=max_gts[i];
        std::unique_ptr<GermlineDiploidIndelLocusInfo> ii(new GermlineDiploidIndelLocusInfo(dopt.gvcf, sampleCount, indelKey,indelData,dindel,indelReportInfo));
        overlap.process(std::move(ii));
    }

    overlap.flush();

    BOOST_REQUIRE_EQUAL(next->the_indels.size(),4u);
}


BOOST_AUTO_TEST_SUITE_END()
