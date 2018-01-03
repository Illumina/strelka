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

#include <iosfwd>

#include <cmath>
#include <cstdint>


/// organizes variant mapq metrics
///
struct MapqTracker
{
    void
    add(const uint8_t mapq)
    {
        count++;
        sumSquare += (mapq*mapq);
        if (mapq==0) zeroCount++;
    }

    void
    clear()
    {
        count = 0;
        zeroCount = 0;
        sumSquare = 0;
    }

    void
    merge(
        const MapqTracker& rhs)
    {
        count += rhs.count;
        zeroCount += rhs.zeroCount;
        sumSquare += rhs.sumSquare;
    }

    double
    getRMS() const
    {
        if (count==0) return 0.;
        return std::sqrt(sumSquare/count);
    }

    double
    getZeroFrac() const
    {
        if (count==0) return 0;
        return (zeroCount/static_cast<double>(count));
    }

    unsigned count = 0;
    unsigned zeroCount = 0;
    double sumSquare = 0;
};

std::ostream&
operator<<(std::ostream& os, const MapqTracker& mapq);
