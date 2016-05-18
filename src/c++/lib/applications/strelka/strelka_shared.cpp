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

///
/// \author Chris Saunders
///

#include "position_somatic_snv_strand_grid.hh"
#include "somatic_indel_grid.hh"
#include "strelka_shared.hh"
#include "strelkaScoringFeatures.hh"
#include "blt_util/blt_exception.hh"

#include <sstream>



strelka_deriv_options::
strelka_deriv_options(
    const strelka_options& opt,
    const reference_contig_segment& ref)
    : base_t(opt,ref)
    , _sscaller_strand_grid(new somatic_snv_caller_strand_grid(opt))
    , _sicaller_grid(new somatic_indel_caller_grid(opt))
{
    if (opt.sfilter.is_depth_filter())
    {
        parse_chrom_depth(opt.sfilter.chrom_depth_file, sfilter.chrom_depth);

        //TODO, verify that chroms match bam chroms
        const std::string& chrom_name(opt.bam_seq_name);
        cdmap_t::const_iterator cdi(sfilter.chrom_depth.find(std::string(chrom_name)));
        if (cdi == sfilter.chrom_depth.end())
        {
            std::ostringstream oss;
            oss << "ERROR: Can't find chromosome: '" << chrom_name << "' in chrom depth file: " << opt.sfilter.chrom_depth_file << "\n";
            throw blt_exception(oss.str().c_str());
        }
        sfilter.expected_chrom_depth=(cdi->second);
        sfilter.max_chrom_depth=(cdi->second*opt.sfilter.max_depth_factor);
        assert(sfilter.max_chrom_depth>=0.);
    }

    sfilter.indelRegionStage=(addPostCallStage(opt.sfilter.indelRegionFlankSize));

    if (opt.isUseSomaticSNVScoring())
    {
        somaticSnvScoringModel.reset(
            new VariantScoringModelServer(
                STRELKA_SNV_SCORING_FEATURES::getFeatureMap(),
                opt.somatic_snv_scoring_model_filename,
                SCORING_CALL_TYPE::SOMATIC,
                SCORING_VARIANT_TYPE::SNV)
        );
    }
    if (opt.isUseSomaticIndelScoring())
    {
        somaticIndelScoringModel.reset(
            new VariantScoringModelServer(
                STRELKA_INDEL_SCORING_FEATURES::getFeatureMap(),
                opt.somatic_indel_scoring_model_filename,
                SCORING_CALL_TYPE::SOMATIC,
                SCORING_VARIANT_TYPE::INDEL)
        );
    }
}



/// dtor required to be in the cpp so that unique ptr can access complete data type
strelka_deriv_options::
~strelka_deriv_options() {}
