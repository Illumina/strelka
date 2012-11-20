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
#ifndef __BAM_DUMPER_HH
#define __BAM_DUMPER_HH

extern "C" {
#include "sam.h"
}


struct bam_dumper {

    bam_dumper(const char* filename,
               const bam_header_t* header);

    ~bam_dumper() {
        if(NULL != _bfp) samclose(_bfp);
    }

    void
    put_record(const bam1_t* brec) {
        samwrite(_bfp,brec);
    }

private:
    samfile_t* _bfp;
};


#endif
