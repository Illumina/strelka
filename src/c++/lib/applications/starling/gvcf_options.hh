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

#pragma once

#include "blt_util/chrom_depth_map.hh"
#include "calibration/featuresetUtil.hh"

#include <string>
#include <vector>


struct gvcf_options
{
    bool
    is_depth_filter() const
    {
        return (!chrom_depth_file.empty());
    }

    bool
    is_gvcf_output() const
    {
        return (not outputPrefix.empty());
    }

    bool
    is_max_ref_rep() const
    {
        return (max_ref_rep >= 0);
    }

    // admin/other:
    std::string chrom_depth_file;
    bool is_skip_header = false;

    // filters:
    double max_depth_factor = 3.;
    bool is_min_gqx = true;
    double min_gqx = 30.;
    bool is_min_homref_gqx = true;
    double min_homref_gqx = 30;
    bool is_max_base_filt = true;
    double max_base_filt = .4; //DPFratio
    bool is_max_snv_sb = true;
    double max_snv_sb = 10;
    bool is_max_snv_hpol = true;
    int max_snv_hpol = -1;
    int max_ref_rep = -1;

    /// LowDepth filter is added if AD sum or DP/DPI is below this value
    unsigned minPassedCallDepth = 3u;

    // blocking scheme:
    //
    // NOTE: bool inserted in this line causes some type of struct initialization error in gcc (4.1 and 4.7)... in blt_options and starling options
    //       as if the struct was sliced. This smells like a compiler error but those should be rare... possible other struct layout error here?
    //
    std::string block_label_prefix = "BLOCKAVG_min";
    unsigned block_percent_tol = 30;
    unsigned block_abs_tol = 3;
    bool is_block_compression = true;

    double block_max_nonref = .2; // what percentage of non-ref bases can a site have and still be included in a non-variant block

    std::string outputPrefix;

    /// file specifying regions that are not compressed in the gvcf:
    std::string nocompress_region_bedfile;

    /// conditional header values to include
    std::vector<std::string> include_headers;
};



struct gvcf_deriv_options
{
    gvcf_deriv_options(
        const gvcf_options& opt,
        const bool isRNA);

    bool
    is_max_depth() const
    {
        return (! chrom_depth.empty());
    }

    std::string block_label;
    cdmap_t chrom_depth;

    const FeatureSet& snvFeatureSet;
    const FeatureSet& snvDevelopmentFeatureSet;
    const FeatureSet& indelFeatureSet;
    const FeatureSet& indelDevelopmentFeatureSet;
};
