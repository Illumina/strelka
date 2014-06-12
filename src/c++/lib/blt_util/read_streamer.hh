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
