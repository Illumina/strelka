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
#ifndef __BAM_STREAMER_HH
#define __BAM_STREAMER_HH


#include "blt_util/bam_record.hh"
//#include "blt_util/read_streamer.hh"

#include <string>



struct bam_streamer {

    explicit
    bam_streamer(const char* filename,
                 const char* region = NULL);
    
    ~bam_streamer();

    bool next();

    const bam_record* get_record_ptr() const {
        if(_is_record_set) return &_brec;
        else               return NULL;
    }

    const char* name() const { return _stream_name.c_str(); }

    unsigned record_no() const { return _record_no; }

    void report_state(std::ostream& os) const;

    const char*
    target_id_to_name(const int32_t tid) const;

    int32_t
    target_name_to_id(const char* seq_name) const;

    const bam_header_t*
    get_header() const { return _bfp->header; }

private:
    bool _is_record_set;
    unsigned _record_no;
    std::string _stream_name;

    bool _is_region;
    std::string _region;

    samfile_t* _bfp;
    bam_index_t* _bidx;
    bam_iter_t _biter;

    bam_record _brec;
};


#endif
