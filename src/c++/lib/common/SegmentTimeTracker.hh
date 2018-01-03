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

///
/// \author Chris Saunders
///

#pragma once

#include "blt_util/time_util.hh"

#include "boost/utility.hpp"


#include <iosfwd>



/// simple segment tracker and reporter
struct SegmentTimeTracker : private boost::noncopyable
{
    explicit
    SegmentTimeTracker(const std::string& outputFile);

    ~SegmentTimeTracker();

    void
    start()
    {
        segmentTime.clear();
        segmentTime.resume();
    }

    void
    stop();

    CpuTimes
    getLastSegmentTime() const
    {
        return segmentTime.getTimes();
    }

private:
    std::ostream* _osPtr;
    TimeTracker segmentTime;
};
