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

/// \author Chris Saunders
///

#include "starling_streams.hh"

#include <iostream>



starling_streams::
starling_streams(const starling_options& opt,
                 const prog_info& pinfo,
                 const bam_header_t* const header,
                 const sample_info& ssi)
    : base_t(opt,pinfo,ssi)
{
    if (opt.gvcf.is_gvcf_output())
    {
        _gvcf_osptr[0] = initialize_gvcf_file(opt,pinfo,opt.gvcf.out_file,header,_gvcf_osptr_auto[0]);
    }

    if (opt.is_realigned_read_file)
    {
        _realign_bam_ptr[0].reset(initialize_realign_bam(opt.is_clobber,pinfo,opt.realigned_read_filename,"realigned-read BAM",header));
    }
}
