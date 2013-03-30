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

#pragma once

#include "blt_util/blt_types.hh"


/// \brief base for objects designed to perform work in a single pass over a position range
///
/// Work progress is communicated via the process_pos() method. This base class is designed to
/// link the worker object with the stage_manager object
///
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
