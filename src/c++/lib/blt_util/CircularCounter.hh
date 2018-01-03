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

/// \file
/// \author Chris Saunders
///

#pragma once

#include <cassert>

#include <vector>


/// A circular buffer of fixed size, S
///
/// - true/false values can be pushed in
/// - total true count among the last S pushes can be queried at any point
///    - count() is O(1) operation
///
struct CircularCounter
{
    CircularCounter(
        const unsigned initSize) :
        _count(0),
        _headPos(0),
        _data(initSize,false)
    {}

    void
    push(const bool val)
    {
        if (_data[_headPos])
        {
            if (!val)
            {
                assert(_count>0);
                _count--;
            }
        }
        else
        {
            if (val) _count++;
        }
        _data[_headPos] = val;
        _headPos = nextPos();
    }

    unsigned
    count() const
    {
        return _count;
    }

    unsigned
    size() const
    {
        return _data.size();
    }

private:

    unsigned
    nextPos() const
    {
        const unsigned pos(_headPos+1);
        if (pos>=size()) return 0;
        return pos;
    }

    unsigned _count;
    unsigned _headPos;
    std::vector<bool> _data;
};
