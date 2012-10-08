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
///
/// \author Chris Saunders
///

#ifndef __POS_PROCESSOR_BASE_H
#define __POS_PROCESSOR_BASE_H

#include "blt_util/blt_types.hh"

//
//
//

struct pos_processor_base {

    pos_processor_base() 
        : _is_skip_process_pos(false) {}

    virtual
    ~pos_processor_base() {}

    void
    check_process_pos(const int stage_no,
                      const pos_t pos) { 
        if(_is_skip_process_pos) return;
        process_pos(stage_no,pos);
    }

    virtual
    void
    process_pos(const int stage_no,
                const pos_t pos) = 0;

protected:
    mutable bool _is_skip_process_pos;
};


#endif
