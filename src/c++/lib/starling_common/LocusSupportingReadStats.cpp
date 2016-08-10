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

#include "LocusSupportingReadStats.hh"

#include <iostream>



std::ostream&
operator<<(
    std::ostream& os,
    const LocusSupportingReadStats& lsrs)
{
    const unsigned ac(lsrs.getAltCount());
    os << "altCount: " << ac << "\n";

    const auto& fwd(lsrs.getCounts(true));
    os << "FWD: ";
    for (unsigned alleleIndex(0); alleleIndex<ac; ++alleleIndex)
    {
        os << fwd.confidentAlleleCount(alleleIndex) << ",";
    }
    os << fwd.nonConfidentCount << "\n";

    const auto& rev(lsrs.getCounts(false));
    os << "REV: ";
    for (unsigned alleleIndex(0); alleleIndex<ac; ++alleleIndex)
    {
        os << rev.confidentAlleleCount(alleleIndex) << ",";
    }
    os << rev.nonConfidentCount << "\n";

    os << "depth: " << lsrs.tier1Depth << "\n";
    return os;
}
