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
#ifndef __VCF_STREAMER_HH
#define __VCF_STREAMER_HH

#include "blt_util/bam_util.hh"
#include "blt_util/tabix_util.hh"
#include "blt_util/vcf_record.hh"


struct vcf_streamer {

    // optionally provide a BAM header to validatr vcf chromosome names against
    //
    explicit
    vcf_streamer(const char* filename,
                 const char* region = NULL,
                 const bam_header_t* bh = NULL);
    
    ~vcf_streamer();

    bool next(const bool is_indel_only=false);

    const vcf_record* get_record_ptr() const {
        if(_is_record_set) return &_vcfrec;
        else               return NULL;
    }

    const char* name() const { return _stream_name.c_str(); }

    unsigned record_no() const { return _record_no; }

    void report_state(std::ostream& os) const;

    //const bam_header_t*
    //get_header() const { return _bfp->header; }

private:
    bool _is_record_set;
    bool _is_stream_end;
    unsigned _record_no;
    std::string _stream_name;

    tabix_t* _tfp;
    ti_iter_t _titer;

    vcf_record _vcfrec;
};


#endif
