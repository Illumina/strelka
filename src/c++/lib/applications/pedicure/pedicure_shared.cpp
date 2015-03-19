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

#include "inovo_shared.hh"
#include "blt_util/blt_exception.hh"
#include "blt_util/chrom_depth_map.hh"

#include <cassert>

#include <sstream>


inovo_deriv_options::
inovo_deriv_options(
    const inovo_options& opt,
    const reference_contig_segment& ref)
    : base_t(opt,ref)
{
    if (opt.dfilter.is_depth_filter())
    {
        parse_chrom_depth(opt.dfilter.chrom_depth_file, dfilter.chrom_depth);

        //TODO, verify that chroms match bam chroms
        const std::string& chrom_name(opt.bam_seq_name);
        cdmap_t::const_iterator cdi(dfilter.chrom_depth.find(std::string(chrom_name)));
        if (cdi == dfilter.chrom_depth.end())
        {
            std::ostringstream oss;
            oss << "ERROR: Can't find chromosome: '" << chrom_name << "' in chrom depth file: " << opt.dfilter.chrom_depth_file << "\n";
            throw blt_exception(oss.str().c_str());
        }
        dfilter.max_depth=(cdi->second*opt.dfilter.max_depth_factor);
        assert(dfilter.max_depth>=0.);
    }
}

/// dtor required to be in the cpp so that unique ptr can access complete data type
inovo_deriv_options::
~inovo_deriv_options() {}
