// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
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

#include "blt_util/bam_dumper.hh"
#include "blt_util/log.hh"

#include <cstdlib>

#include <iostream>



bam_dumper::
bam_dumper(const char* filename,
           const bam_header_t* header) {

    _bfp = samopen(filename, "wb", header);

    if(NULL == _bfp) {
        log_os << "ERROR: Failed to open output BAM file: " << filename << "\n";
        exit(EXIT_FAILURE);
    }
}
