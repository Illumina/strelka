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

#include <cassert>
#include <cstdint>


namespace STAR_DIINDEL
{
enum index_t
{
    NOINDEL,
    HOM,
    HET,
    SIZE
};

inline
const char*
label(const unsigned idx)
{
    switch (idx)
    {
    case NOINDEL:
        return "ref";
    case HOM:
        return "hom";
    case HET:
        return "het";
    default:
        return "xxx";
    }
}

inline
const char*
get_gt_label(const unsigned idx)
{
    switch (idx)
    {
    case NOINDEL:
        return "0/0";
    case HOM:
        return "1/1";
    case HET:
        return "0/1";
    default:
        assert(false && "Unknown Indel GT");
        return nullptr;
    }
}

inline
uint8_t
get_allele(
    const unsigned idx,
    const unsigned chromidx)
{
    assert(idx<SIZE);
    assert(chromidx<2);
    static const uint8_t res[SIZE][2] = {{0,0},{1,1},{0,1}};
    return res[idx][chromidx];
}

}
