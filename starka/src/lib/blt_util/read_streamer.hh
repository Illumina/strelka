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
#ifndef __READ_STREAMER_HH
#define __READ_STREAMER_HH

#include "blt_util/read_record.hh"

#include <boost/utility.hpp>

#include <iosfwd>


// unified interface between bam and export readers:
//
struct read_streamer : private boost::noncopyable {

    virtual ~read_stream() {}

    virtual bool next() = 0;
    //    virtual const read_record* get_record_ptr() const = 0;
    virtual const char* name() const = 0;
    virtual unsigned record_no() const = 0;
    virtual void report_state(std::ostream& os) const = 0;
};



#endif
