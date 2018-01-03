//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

/// \file

/// \author Chris Saunders
///

#pragma once

#include "blt_util/read_record.hh"

#include <boost/utility.hpp>

#include <iosfwd>


// unified interface between bam and export readers:
//
struct read_streamer : private boost::noncopyable
{
    virtual ~read_streamer() {}

    virtual bool next() = 0;
    //    virtual const read_record* get_record_ptr() const = 0;
    virtual const char* name() const = 0;
    virtual unsigned record_no() const = 0;
    virtual void report_state(std::ostream& os) const = 0;
};
