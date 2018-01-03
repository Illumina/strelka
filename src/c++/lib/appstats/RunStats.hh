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

#include "boost/serialization/nvp.hpp"
#include "boost/serialization/vector.hpp"

#include <cassert>
#include <cstdint>

#include <iosfwd>
#include <vector>


struct RunStatsData
{
    RunStatsData() = default;

    void
    merge(const RunStatsData& rhs)
    {
        lifeTime.merge(rhs.lifeTime);
        candidateIndels += rhs.candidateIndels;
        nonCandidateIndels += rhs.nonCandidateIndels;
    }

    void
    report(std::ostream& os) const;

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& BOOST_SERIALIZATION_NVP(lifeTime);
        ar& BOOST_SERIALIZATION_NVP(candidateIndels);
        ar& BOOST_SERIALIZATION_NVP(nonCandidateIndels);
    }

    /// Total wall-time of each (single-thread) process, summed together
    CpuTimes lifeTime;

    /// Total indels reaching candidate status in the report range and (if defined) call regions
    unsigned long candidateIndels = 0;

    /// Total indels failing to reach candidate status in the report range and (if defined) call regions
    unsigned long nonCandidateIndels = 0;
};

BOOST_CLASS_IMPLEMENTATION(RunStatsData, boost::serialization::object_serializable)

struct RunStats
{
    void
    load(const char* filename);

    void
    save(std::ostream& os) const;

    void
    save(const char* filename) const;

    void
    report(const char* filename) const;

    void
    merge(const RunStats& rhs)
    {
        runStatsData.merge(rhs.runStatsData);
    }

    RunStatsData runStatsData;
};

BOOST_CLASS_IMPLEMENTATION(RunStats, boost::serialization::object_serializable)
