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

#include "VariantPhaser.hh"

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

BOOST_AUTO_TEST_SUITE( variantPhaserTest )



BOOST_AUTO_TEST_CASE( simplePhasingTest )
{
//    reference_contig_segment rcs;
//    rcs.seq() = "ACGGGGTTGGACGATGCTACGATCGATCGCGTACCTACGATCGACTACGACTGCGACGATCGACGATCGACGATCGATCGATCGACGTACGACACGTACGATCGATCGATCGATCGACTCGATCAGCTCATGCATCG";
//
//    starling_options opt;
//    opt.is_user_genome_size = true;
//    opt.user_genome_size = rcs.seq().size();
//
//    starling_deriv_options dopt(opt);
//
//    std::shared_ptr<variant_pipe_stage_base> next(new dummy_variant_sink);
//
//    const unsigned sampleCount(1);
//    VariantPhaser phaser(opt, sampleCount, next);
//
//    IndelKey indelKey;
//    const IndelData indelData(1,indelKey);
//
//    indelKey.pos=6;
//    indelKey.type=INDEL::INDEL;
//    indelKey.deletionLength=2;
//
//    indelData.getSampleData(0).haplotypeId = 0;
//    std::unique_ptr<GermlineDiploidIndelLocusInfo> ii(new GermlineDiploidIndelLocusInfo(dopt.gvcf, sampleCount, -1));
//    ii->addAltIndelAllele(indelKey, indelData);
//    phaser.process(std::move(ii));
//    phaser.flush();
}

BOOST_AUTO_TEST_SUITE_END()
