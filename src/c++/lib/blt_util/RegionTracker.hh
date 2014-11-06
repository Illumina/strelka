// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#pragma once

#include "blt_util/known_pos_range2.hh"

#include <iosfwd>
#include <set>


/// sort pos range using end_pos as the primary sort key
struct PosRangeEndSort
{
    bool
    operator()(
        const known_pos_range2& lhs,
        const known_pos_range2& rhs) const
    {
        if (lhs.end_pos() < rhs.end_pos()) return true;
        if (lhs.end_pos() == rhs.end_pos())
        {
            if (lhs.begin_pos() < rhs.begin_pos()) return true;
        }
        return false;
    }
};


/// facilitate 'rolling' region tracking and position intersect queries
///
struct RegionTracker
{
    bool
    isInRegion(const unsigned pos) const;

    /// add region
    ///
    /// any overlaps with existing regions in the tracker will be collapsed
    void
    addRegion(known_pos_range2 range);

    /// remove all regions which end (inclusive) before pos+1
    void
    removeToPos(const unsigned pos);

    // debug util
    void
    dump(std::ostream& os) const;

private:
    std::set<known_pos_range2,PosRangeEndSort> _regions;
};
