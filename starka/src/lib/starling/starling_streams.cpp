// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright 2009 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
//

/// \file

/// \author Chris Saunders
///

#include "starling_streams.hh"

#include <iostream>


const sample_info starling_sample_info;



starling_streams::
starling_streams(const starling_options& opt,
                 const prog_info& pinfo,
                 const bam_header_t* const header)
    : base_t(opt,pinfo,starling_sample_info) {

    if(opt.is_bindel_diploid_file){
        _bindel_diploid_osptr[0].reset(initialize_bindel_file(opt,pinfo,opt.bindel_diploid_filename));
    }

    if(opt.is_realigned_read_file) {
        _realign_bam_ptr[0].reset(initialize_realign_bam(opt.is_clobber,pinfo,opt.realigned_read_filename,"realigned-read BAM",header));
    }
}
