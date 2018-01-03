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

#include <cstdint>


struct SiteNoise
{
    void
    clear()
    {
        total = 0;
        noise = 0;
        noise2 = 0;
    }

    double
    nfrac() const
    {
        return static_cast<double>(noise)/total;
    }

    double
    n2frac() const
    {
        return static_cast<double>(noise2)/total;
    }

    uint16_t total = 0;
    uint16_t noise = 0;
    uint16_t noise2 = 0;
};


void
set_noise_from_vcf(
    const char* line,
    SiteNoise& sn);

