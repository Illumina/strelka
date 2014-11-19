// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#pragma once

#include "blt_util/chrom_depth_map.hh"

#include <string>


struct gvcf_options
{
    bool
    is_depth_filter() const
    {
        return (! chrom_depth_file.empty());
    }

    bool
    is_gvcf_output() const
    {
        return (! out_file.empty());
    }


    // admin/other:
    std::string chrom_depth_file;
    bool is_skip_header = false;

    // filters:
    double max_depth_factor = 3.;
    bool is_min_gqx = true;
    double min_gqx = 30.;
    bool is_max_base_filt = true;
    double max_base_filt = .4; //DPFratio
    bool is_max_snv_sb = true;
    double max_snv_sb = 10;
    bool is_max_snv_hpol = true;
    int max_snv_hpol = -1;
    bool is_max_ref_rep = true;
    int max_ref_rep = -1;

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
    std::string out_file;

    /// file specifying regions that are not compressed in the gvcf:
    std::string nocompress_region_bedfile;
};



struct gvcf_deriv_options
{
    gvcf_deriv_options(
        const gvcf_options& opt,
        const std::string& bam_seq_name);

    bool
    is_max_depth() const
    {
        return (! chrom_depth.empty());
    }

    std::string block_label;
    double max_depth = 0;
    cdmap_t chrom_depth;
};
