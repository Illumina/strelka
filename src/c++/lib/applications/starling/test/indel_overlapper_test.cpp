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

    void process(std::unique_ptr<GermlineSiteCallInfo> si) override
    {
        the_sites.push_back(std::move(si));
    }
    void process(std::unique_ptr<GermlineIndelCallInfo> ii) override
    {
        the_indels.push_back(std::move(ii));
    }

    std::vector<std::unique_ptr<GermlineSiteCallInfo>> the_sites;
    std::vector<std::unique_ptr<GermlineIndelCallInfo>> the_indels;
};



BOOST_AUTO_TEST_SUITE( indel_overlapper_test )

/// add a single indel to the indel_overlapper, this is
/// a simple demo to unit test the overlapper
///
BOOST_AUTO_TEST_CASE( simple_indel_test )
{
    // fake various high-level data structures with as many defaults as possible
    starling_options opt;

    reference_contig_segment rcs;
    rcs.seq() = "ACGGGGTTGGACGATGCTACGATCGATCGCGTACCTACGATCGACTACGACTGCGACGATCGACGATCGACGATCGATCGATCGACGTACGACACGTACGATCGATCGATCGATCGACTCGATCAGCTCATGCATCG";

    gvcf_deriv_options gvcf_dopt(opt.gvcf, "chr1");

    ScoringModelManager cm(opt,gvcf_dopt);

    std::shared_ptr<variant_pipe_stage_base> next(new dummy_variant_sink);
    indel_overlapper overlap(cm, rcs, next);

    indel_key ik;
    const IndelData id(1,ik);
    const GermlineDiploidIndelSimpleGenotypeInfoCore dindel;
    const starling_indel_report_info iri;
    const starling_indel_sample_report_info isri;

    ik.pos=6;
    ik.type=INDEL::DELETE;
    ik.length=2;

    std::unique_ptr<GermlineDiploidIndelCallInfo> ii(new GermlineDiploidIndelCallInfo(ik,id,dindel,iri,isri));
    overlap.process(std::move(ii));

    overlap.flush();
    // just testing that we don't assert for now...
}


/// Reduced recreation of STARKA-393 failure:
BOOST_AUTO_TEST_CASE( conflicting_indel_test )
{
    // fake various high-level data structures with as many defaults as possible
    starling_options opt;

    reference_contig_segment rcs;
    rcs.seq() = "ACGGGGTTGGACGATGCTACGATCGATCGCGTACCTACGATCGACTACGACTGCGACGATCGACGATCGACGATCGATCGATCGACGTACGACACGTACGATCGATCGATCGATCGACTCGATCAGCTCATGCATCG";

    gvcf_deriv_options gvcf_dopt(opt.gvcf, "chr1");

    ScoringModelManager cm(opt,gvcf_dopt);

    std::shared_ptr<variant_pipe_stage_base> next(new dummy_variant_sink);
    indel_overlapper overlap(cm, rcs, next);

    indel_key iks[] =
    {
        indel_key(10,INDEL::DELETE,10),
        indel_key(15,INDEL::DELETE,30),
        indel_key(20,INDEL::DELETE,1),
        indel_key(25,INDEL::DELETE,1),
    };

    int max_gts[] = { 2,0,2,2 };

    indel_key ik;
    const IndelData id(1,ik);
    GermlineDiploidIndelSimpleGenotypeInfoCore dindel;
    dindel.is_forced_output = true;
    const starling_indel_report_info iri;
    const starling_indel_sample_report_info isri;

    const unsigned keyCount(sizeof(iks)/sizeof(indel_key));
    for (unsigned i(0); i<keyCount; ++i)
    {
        ik=iks[i];
        dindel.max_gt=max_gts[i];
        dindel.max_gt_poly=max_gts[i];
        std::unique_ptr<GermlineDiploidIndelCallInfo> ii(new GermlineDiploidIndelCallInfo(ik,id,dindel,iri,isri));
        overlap.process(std::move(ii));
    }

    overlap.flush();
    // just testing that we don't assert for now...
}


/// Additional test added in light of STARKA-393 failure:
BOOST_AUTO_TEST_CASE( conflicting_indel_test2 )
{
    // fake various high-level data structures with as many defaults as possible
    starling_options opt;

    reference_contig_segment rcs;
    rcs.seq() = "ACGGGGTTGGACGATGCTACGATCGATCGCGTACCTACGATCGACTACGACTGCGACGATCGACGATCGACGATCGATCGATCGACGTACGACACGTACGATCGATCGATCGATCGACTCGATCAGCTCATGCATCG";

    gvcf_deriv_options gvcf_dopt(opt.gvcf, "chr1");

    ScoringModelManager cm(opt,gvcf_dopt);

    std::shared_ptr<dummy_variant_sink> next(new dummy_variant_sink);
    indel_overlapper overlap(cm, rcs, next);

    indel_key iks[] =
    {
        indel_key(10,INDEL::DELETE,10),
        indel_key(12,INDEL::DELETE,1),
        indel_key(15,INDEL::DELETE,1),
        indel_key(18,INDEL::DELETE,1),
    };

    int max_gts[] = { 2,2,0,2 };

    indel_key ik;
    const IndelData id(1,ik);
    GermlineDiploidIndelSimpleGenotypeInfoCore dindel;
    dindel.is_forced_output = true;
    const starling_indel_report_info iri;
    const starling_indel_sample_report_info isri;

    const unsigned keyCount(sizeof(iks)/sizeof(indel_key));
    for (unsigned i(0); i<keyCount; ++i)
    {
        ik=iks[i];
        dindel.max_gt=max_gts[i];
        dindel.max_gt_poly=max_gts[i];
        std::unique_ptr<GermlineDiploidIndelCallInfo> ii(new GermlineDiploidIndelCallInfo(ik,id,dindel,iri,isri));
        overlap.process(std::move(ii));
    }

    overlap.flush();

    BOOST_REQUIRE_EQUAL(next->the_indels.size(),4u);
}


BOOST_AUTO_TEST_SUITE_END()
