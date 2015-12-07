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

#pragma once

#include "blt_util/blt_types.hh"
#include "blt_util/RangeMap.hh"

#include <cassert>



/// simple map of position to depth
///
/// assumes that a narrow list of positions is maintained so that
/// array based lookup optimizations can be used
///
struct depth_buffer
{
    unsigned
    val(const pos_t pos) const
    {
        return _data.getConstRefDefault(pos,0);
    }

    void
    inc(const pos_t pos)
    {
        _data.getRef(pos) += 1;
    }

    void
    clear_pos(const pos_t pos)
    {
        if (_data.isKeyPresent(pos)) _data.erase(pos);
    }

    /// return true if buffered depth exceeds depth in [begin,end]
    bool
    is_range_ge_than(const pos_t begin,
                     const pos_t end,
                     const unsigned depth) const
    {
        assert(begin <= end);
        for (pos_t i(begin); i<=end; ++i)
        {
            if (val(i) >= depth) return true;
        }
        return false;
    }

private:
    typedef RangeMap<pos_t,unsigned> count_t;
    count_t _data;
};
