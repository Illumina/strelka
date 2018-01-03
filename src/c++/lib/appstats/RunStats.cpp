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

#include "RunStats.hh"

#include "boost/archive/xml_iarchive.hpp"
#include "boost/archive/xml_oarchive.hpp"

#include <fstream>
#include <iostream>



void
RunStatsData::
report(std::ostream& os) const
{
    os << "TotalHours\t";
    lifeTime.reportHr(os);
    os << "\n";
    os << "\n";
    os << "CallRegionCandidateIndels\t" << candidateIndels << "\n";
    os << "CallRegionNonCandidateIndels\t" << nonCandidateIndels << "\n";
}



void
RunStats::
load(const char* filename)
{
    assert(nullptr != filename);
    std::ifstream ifs(filename);
    boost::archive::xml_iarchive ia(ifs);
    ia >> BOOST_SERIALIZATION_NVP(runStatsData);
}



void
RunStats::
save(std::ostream& os) const
{
    boost::archive::xml_oarchive oa(os);
    oa << BOOST_SERIALIZATION_NVP(runStatsData);
}



void
RunStats::
save(const char* filename) const
{
    assert(nullptr != filename);
    std::ofstream ofs(filename);
    save(ofs);
}



void
RunStats::
report(const char* filename) const
{
    assert(nullptr != filename);
    std::ofstream ofs(filename);
    ofs << "StrelkaRunStatsReport\n";
    runStatsData.report(ofs);
}
