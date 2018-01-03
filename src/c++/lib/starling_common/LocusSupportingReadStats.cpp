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

#include "LocusSupportingReadStats.hh"

#include <iostream>



std::ostream&
operator<<(
    std::ostream& os,
    const LocusSupportingReadStats& lsrs)
{
    const unsigned altAlleleCount(lsrs.getAltCount());
    const unsigned fullAlleleCount(altAlleleCount+1);
    os << "altAlleleCount: " << altAlleleCount << "\n";

    auto dumpStrandCounts = [&](const bool isFwdStrand)
    {
        const auto& counts(lsrs.getCounts(isFwdStrand));
        os << "AD: ";
        for (unsigned alleleIndex(0); alleleIndex < fullAlleleCount; ++alleleIndex)
        {
            if (alleleIndex>0) os << ',';
            os << counts.confidentAlleleCount(alleleIndex);
        }
        os << " uncertain: " << counts.nonConfidentCount << "\n";
    };

    os << "FwdStrand: ";
    dumpStrandCounts(true);

    os << "RevStrand: ";
    dumpStrandCounts(false);

    return os;
}
