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


namespace MONOGT
{
enum index_t
{
    A,
    C,
    G,
    T,
    SIZE
};

inline
const char*
label(const unsigned idx)
{
    switch (idx)
    {
    case A:
        return "A";
    case C:
        return "C";
    case G:
        return "G";
    case T:
        return "T";
    default:
        return "X";
    }
}

// the lhood function is no longer so general that these values can actually be changed...
inline
double
expect(const int base_id,
       const int gt)
{

    static const unsigned N_BASE(4);

    static const double ex[SIZE][N_BASE] = {{ 1.0, 0.0, 0.0, 0.0},
        { 0.0, 1.0, 0.0, 0.0},
        { 0.0, 0.0, 1.0, 0.0},
        { 0.0, 0.0, 0.0, 1.0}
    };

    return ex[gt][base_id];
}
}
