// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
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

#include "ploidy_util.hh"
#include "common/Exceptions.hh"

#include <sstream>



boost::optional<int>
parsePloidyFromBed(const char* line)
{
    boost::optional<int> result;

    if (line == nullptr) return result;

    unsigned tabcount(0);
    while (true)
    {
        if (*line=='\0' || *line=='\n') return result;
        if (*line=='\t') tabcount++;
        line++;
        if (tabcount>=4) break;
    }

    const char* s(line);
    int val = illumina::blt_util::parse_int(s);
    if (s != line) result.reset(val);

    return result;
}



int
parsePloidyFromBedStrict(const char* line)
{
    using namespace illumina::common;
    const auto ploidy = parsePloidyFromBed(line);
    if (! ploidy)
    {
        std::ostringstream oss;
        oss << "ERROR: can't parse ploidy (column 5) from bed record: '" << line << "'\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }
    if ((*ploidy < 0) || (*ploidy > 2))
    {
        std::ostringstream oss;
        oss << "ERROR: parsed unsupported ploidy value (" << *ploidy << ") from bed record: '" << line << "'\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }
    return *ploidy;
}
