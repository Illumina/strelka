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

#include "blt_util/window_util.hh"


/// \brief Manage all data structures which provide statistics averaged over a range of genome positions
///
/// These averaged statistics are usually used for ranges surrounding a variant call to ascertain properties of
/// the region surrounding the variant.
///
struct LocalRegionStats
{
    explicit
    LocalRegionStats(const unsigned regionSize)
        : regionUsedBasecallCount(regionSize)
        , regionUnusedBasecallCount(regionSize)
        , regionSpanningDeletionReadCount(regionSize)
        , regionSubmappedReadCount(regionSize)
    {}

    void
    reset()
    {
        regionUsedBasecallCount.reset();
        regionUnusedBasecallCount.reset();
        regionSpanningDeletionReadCount.reset();
        regionSubmappedReadCount.reset();
    }

    window_average regionUsedBasecallCount;
    window_average regionUnusedBasecallCount;
    window_average regionSpanningDeletionReadCount;
    window_average regionSubmappedReadCount;
};
