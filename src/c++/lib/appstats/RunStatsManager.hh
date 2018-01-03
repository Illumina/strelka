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

#pragma once

#include "RunStats.hh"
#include "blt_util/time_util.hh"

#include "boost/utility.hpp"

#include <iosfwd>
#include <string>


/// \brief Handles all messy real world interaction for the stats module, while the stats module itself just
///        accumulates data
///
struct RunStatsManager : private boost::noncopyable
{
    explicit
    RunStatsManager(
        const std::string& outputFile);

    ~RunStatsManager();

    void
    addCallRegionIndel(const bool isCandidate)
    {
        if (isCandidate)
        {
            runStats.runStatsData.candidateIndels++;
        }
        else
        {
            runStats.runStatsData.nonCandidateIndels++;
        }
    }

private:
    std::ostream* _osPtr;

    /// this object tracks its own lifetime from ctor-to-dtor here, this is used to approximate
    /// program lifetime when RunStatsManager is appropriately scoped
    TimeTracker lifeTime;

    /// runStats is the primary stats data store
    RunStats runStats;
};
