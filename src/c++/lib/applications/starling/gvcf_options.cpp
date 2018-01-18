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

#include "gvcf_options.hh"
#include "germlineVariantEmpiricalScoringFeatures.hh"
#include "rnaVariantEmpiricalScoringFeatures.hh"

#include <sstream>
#include <iostream>



gvcf_deriv_options::
gvcf_deriv_options(
    const gvcf_options& opt,
    const bool isRNA)
    : snvFeatureSet(isRNA ? RNA_SNV_SCORING_FEATURES::getInstance() : GERMLINE_SNV_SCORING_FEATURES::getInstance()),
      snvDevelopmentFeatureSet(isRNA ? RNA_SNV_SCORING_DEVELOPMENT_FEATURES::getInstance() : GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::getInstance()),
      indelFeatureSet(isRNA ? RNA_INDEL_SCORING_FEATURES::getInstance() : GERMLINE_INDEL_SCORING_FEATURES::getInstance()),
      indelDevelopmentFeatureSet(isRNA ? RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::getInstance() : GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::getInstance())
{
    {
        std::ostringstream oss;
        oss << opt.block_label_prefix << opt.block_percent_tol << "p" << opt.block_abs_tol << "a";
        block_label.assign(oss.str());
    }

    if (opt.is_depth_filter() && opt.is_gvcf_output())
    {
        parse_chrom_depth(opt.chrom_depth_file, chrom_depth);
    }
}
