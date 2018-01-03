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

#include <cassert>


struct Tier2Options
{
    void
    validate() const
    {
        assert((randomBaseMatchProb >= 0.) && (randomBaseMatchProb <= 1.));
    }

    /// min mapping quality used for tier2 calling
    int minMappingErrorPhredProb = 0;

    /// number of mismatches allowed in the mismatch density filter window at tier2
    int mismatchDensityFilterMaxMismatchCount = 10;

    /// If true, use reads with unmapped mates for tier2 calling
    bool includeSingletonReads = true;

    /// If true, use non proper-pair reads for tier2 calling
    bool includeAnomalousReads = true;

    /// If true, use tier2 randomBaseMatchProb value
    bool isRandomBaseMatchProb = true;

    /// tier2 version of the value in starling_base_options
    double randomBaseMatchProb = 0.25;
};
