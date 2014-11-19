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

/// \file

/// \author Chris Saunders
///

#include "gvcf_options.hh"
#include "blt_util/blt_exception.hh"

#include <cassert>

#include <sstream>
#include <iostream>



gvcf_deriv_options::
gvcf_deriv_options(
    const gvcf_options& opt,
    const std::string& bam_seq_name)
{
    {
        std::ostringstream oss;
        oss << opt.block_label_prefix << opt.block_percent_tol << "p" << opt.block_abs_tol << "a";
        block_label.assign(oss.str());
    }

    if (opt.is_depth_filter() && opt.is_gvcf_output())
    {
        parse_chrom_depth(opt.chrom_depth_file, chrom_depth);

        //TODO, verify that chroms match bam chroms
        const std::string& chrom_name(bam_seq_name);
        cdmap_t::const_iterator cdi(chrom_depth.find(std::string(chrom_name)));
        if (cdi == chrom_depth.end())
        {
            std::ostringstream oss;
            oss << "ERROR: Can't find chromosome: '" << chrom_name << "' in chrom depth file: " << opt.chrom_depth_file << "\n";
            throw blt_exception(oss.str().c_str());
        }
        max_depth=(cdi->second*opt.max_depth_factor);
        norm_depth=(cdi->second);
        assert(max_depth>=0.);
        assert(norm_depth>=0.);
    }
}
