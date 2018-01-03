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

/// variation on the original strawman snv caller -- implements a
/// compile-time specified grid in allele frequency space and requires
/// similar frequency as definition of non-somatic.
///

/// \author Chris Saunders
///

#pragma once

// a simplification of diploid calls down to types relative to the reference:
//
namespace NTYPE
{

enum index_t { REF,
               HOM,
               HET,
               CONFLICT,
               SIZE
             };

inline
const char*
label(const unsigned idx)
{
    switch (idx)
    {
    case REF:
        return "ref";
    case HOM:
        return "hom";
    case HET:
        return "het";
    case CONFLICT:
        return "conflict";
    default:
        return "xxx";
    }
}
}
