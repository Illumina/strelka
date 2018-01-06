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

#include "RunStatsManager.hh"
#include "common/Exceptions.hh"

#include <fstream>
#include <iostream>
#include <sstream>



RunStatsManager::
RunStatsManager(
    const std::string& outputFile)
    : _osPtr(nullptr)
{
    if (outputFile.empty()) return;
    _osPtr = new std::ofstream(outputFile.c_str());
    if (! *_osPtr)
    {
        std::ostringstream oss;
        oss << "Can't open output file: '" << outputFile << "'";
        BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
    }

    lifeTime.resume();
}



RunStatsManager::
~RunStatsManager()
{
    if (_osPtr != nullptr)
    {
        lifeTime.stop();
        runStats.runStatsData.lifeTime=lifeTime.getTimes();
        runStats.save(*_osPtr);
        delete _osPtr;
    }
}
