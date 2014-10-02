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

#include "position_somatic_snv.hh"
#include "position_somatic_snv_strand_grid.hh"
#include "somatic_indel.hh"
#include "somatic_indel_grid.hh"
#include "strelka_shared.hh"



strelka_deriv_options::
strelka_deriv_options(
    const strelka_options& opt,
    const reference_contig_segment& ref)
    : base_t(opt,ref)
    , _sscaller_strand_grid(new somatic_snv_caller_strand_grid(opt,pdcaller()))
    , _sicaller_grid(new somatic_indel_caller_grid(opt,incaller()))
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
        sfilter.max_depth=(cdi->second*opt.sfilter.max_depth_factor);
        assert(sfilter.max_depth>=0.);
    }

    sfilter.indelRegionStage=(addPostCallStage(opt.sfilter.indelRegionFlankSize));
}

/// dtor required to be in the cpp so that unique ptr can access complete data type
strelka_deriv_options::
~strelka_deriv_options() {}
