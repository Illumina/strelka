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

#include "test_config.h"


#include "boost/test/unit_test.hpp"

#include "../ScoringModelManager.hh"



BOOST_AUTO_TEST_SUITE( cmodel )


BOOST_AUTO_TEST_CASE( test_cmodel_gqx )
{
#if 0
    // disable for now. These don't fail in unit test form, only integrated
    starling_options opts;
    opts.calibration_models_filename = TEST_CONFIG_PATH;
    opts.is_clobber = true;
    opts.min_paired_align_score=20;
    opts.min_single_align_score=20;
    opts.min_qscore=17 ;
    opts.bsnp_ssd_no_mismatch=0.35 ;
    opts.bsnp_ssd_one_mismatch=0.6 ;
    opts.min_vexp=0.25;

    opts.calibration_model = "QScoreHPDRE100_v4";

    std::vector<indel_info> indels;
    indels.resize(2);


    gvcf_deriv_options gvcf_options(opts.gvcf, "chr2");
    gvcf_options.chrom_depth["chr2"] = 30.734;
    calibration_models cm(opts, gvcf_options);


    IndelKey ik(233030782, INDEL::INSERT, 4);
    starling_diploid_indel  dindel;
    dindel.ploidy=2;
    dindel.max_gt=STAR_DIINDEL::HET;
    dindel.indel_qphred=219;
    dindel.max_gt_qphred=219;
    dindel.max_gt_poly=2;
    dindel.max_gt_poly_qphred=236;
    dindel.is_forced_output=0;
    dindel.is_zero_coverage=0;
    std::initializer_list<double> init({1.34748e-22,2.59853e-24,1});
    std::copy(init.begin(), init.end(), dindel.pprob);

    starling_indel_report_info iri;
    iri.desc="4I";
    iri.ref_seq="----";
    iri.indel_seq="AATA";
    iri.vcf_ref_seq="C";
    iri.vcf_indel_seq="CAATA";
    iri.ref_upstream="CTGTCTCTAC";
    iri.ref_downstream="AATAAATAAA";
    iri.repeat_unit="AATA";
    iri.ref_repeat_count=11;
    iri.indel_repeat_count=12;
    iri.ihpol=5;
    iri.it=INDEL::INSERT;
    starling_indel_sample_report_info isri;
    isri.n_q30_ref_reads=1;
    isri.n_q30_indel_reads=6;
    isri.n_q30_alt_reads=4;
    isri.n_other_reads=17;
    isri.depth=28;

    indel_info& ii(indels[0]);
    ii.init(ik.pos,ik,dindel,iri,isri);
    ii.imod.max_gt = DIGT::AA;
    cigar_to_apath("1M4I", ii.imod.cigar);
    ii.imod.is_overlap = true;


    cm.classify_indel(ii, ii.imod);
    //BOOST_CHECK_EQUAL(ii.imod.gqx, 219);
    BOOST_CHECK_EQUAL(ii.imod.filters.none(), true);



    IndelKey ik2(233030782, INDEL::INSERT, 8);

    starling_diploid_indel  dindel2;

    dindel2.ploidy=2;
    dindel2.max_gt=STAR_DIINDEL::HET;
    dindel2.indel_qphred=120;
    dindel2.max_gt_qphred=120;
    dindel2.max_gt_poly=2;
    dindel2.max_gt_poly_qphred=163;
    dindel2.is_forced_output=0;
    dindel2.is_zero_coverage=0;
    std::initializer_list<double> init2({1.03625e-12,6.35264e-36,1});
    std::copy(init2.begin(), init2.end(), dindel2.pprob);

    starling_indel_report_info iri2;
    iri2.desc="8I";
    iri2.ref_seq="--------";
    iri2.indel_seq="AATAAATA";
    iri2.vcf_ref_seq="C";
    iri2.vcf_indel_seq="CAATAAATA";
    iri2.ref_upstream="CTGTCTCTAC";
    iri2.ref_downstream="AATAAATAAA";
    iri2.repeat_unit="AATA";
    iri2.ref_repeat_count=11;
    iri2.indel_repeat_count=13;
    iri2.ihpol=5;
    iri2.it=INDEL::INSERT;

    starling_indel_sample_report_info isri2;
    isri2.n_q30_ref_reads=1;
    isri2.n_q30_indel_reads=4;
    isri2.n_q30_alt_reads=6;
    isri2.n_other_reads=17;
    isri2.depth=28;

    indel_info& ii2(indels[1]);
    ii2.init(ik2.pos,ik2,dindel2,iri2,isri2);
    ii2.imod.max_gt = DIGT::AA;
    cigar_to_apath("1M8I", ii2.imod.cigar);
    ii2.imod.is_overlap = true;

    cm.classify_indel(ii2, ii2.imod);
    // BOOST_CHECK_EQUAL(ii2.imod.gqx, 120);
    BOOST_CHECK_EQUAL(ii2.imod.filters.none(), true);

#endif


}


BOOST_AUTO_TEST_SUITE_END()
