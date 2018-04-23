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

#include "VariantOverlapResolver.hh"

#include "ScoringModelManager.hh"


struct DummyVariantSink : public variant_pipe_stage_base
{
    DummyVariantSink() : variant_pipe_stage_base() {}

    void process(std::unique_ptr<GermlineSiteLocusInfo> si) override
    {
        _siteLoci.push_back(std::move(si));
    }
    void process(std::unique_ptr<GermlineIndelLocusInfo> ii) override
    {
        _indelLoci.push_back(std::move(ii));
    }

    std::vector<std::unique_ptr<GermlineSiteLocusInfo>> _siteLoci;
    std::vector<std::unique_ptr<GermlineIndelLocusInfo>> _indelLoci;
};



BOOST_AUTO_TEST_SUITE( indel_overlapper_test_suite )

/// add a single indel to the indel_overlapper, this is
/// a simple demo to unit test the overlapper
///
BOOST_AUTO_TEST_CASE( simple_indel_test )
{
    // fake various high-level data structures with as many defaults as possible
    reference_contig_segment rcs;
    rcs.seq() = "ACGGGGTTGGACGATGCTACGATCGATCGCGTACCTACGATCGACTACGACTGCGACGATCGACGATCGACGATCGATCGATCGACGTACGACACGTACGATCGATCGATCGATCGACTCGATCAGCTCATGCATCG";

    starling_options opt;
    opt.alignFileOpt.alignmentFilenames.push_back("sample.bam");

    const starling_deriv_options dopt(opt);

    ScoringModelManager cm(opt, dopt.gvcf);

    std::shared_ptr<variant_pipe_stage_base> next(new DummyVariantSink);
    VariantOverlapResolver overlap(cm, next);

    IndelKey indelKey;
    const IndelData indelData(1,indelKey);

    indelKey.pos=6;
    indelKey.type=INDEL::INDEL;
    indelKey.deletionLength=2;

    const unsigned sampleCount(1);
    std::unique_ptr<GermlineDiploidIndelLocusInfo> ii(new GermlineDiploidIndelLocusInfo(dopt.gvcf, sampleCount));
    ii->addAltIndelAllele(indelKey, indelData);
    overlap.process(std::move(ii));

    overlap.flush();

    // just testing that we don't assert for now...
    BOOST_REQUIRE(true);
}


BOOST_AUTO_TEST_SUITE_END()
