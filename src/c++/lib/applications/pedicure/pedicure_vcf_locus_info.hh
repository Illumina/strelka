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

#include <cassert>

#include <bitset>
#include <iosfwd>


namespace PEDICURE_VCF_FILTERS
{

enum index_t
{
    // SNVs and indels:
    HighDepth,
    WrongCount,
    // snvs only:
    QDS,
    // indels only:
    QDI,
    Repeat,
    iHpol,
    overlapConflict,
    lowGQX,
    SIZE
};

inline
const char*
get_label(const unsigned idx)
{
    switch (idx)
    {
    case HighDepth:
        return "HighDepth";
    case WrongCount:
        return "WrongCount";
    case Repeat:
        return "Repeat";
    case iHpol:
        return "iHpol";
    case lowGQX:
        return "lowGQX";
    case overlapConflict:
        return "overlapConflict";
    default:
        assert(false && "Unknown vcf filter id");
        return nullptr;
    }
}
}

struct pedicure_shared_modifiers
{
    pedicure_shared_modifiers()
    {
        clear();
    }

    void
    set_filter(const PEDICURE_VCF_FILTERS::index_t i)
    {
        filters.set(i);
    }

    void
    write_filters(std::ostream& os) const;

    void
    clear()
    {
        filters.reset();
    }

    std::bitset<PEDICURE_VCF_FILTERS::SIZE> filters;
};


std::ostream& operator<<(std::ostream& os,const pedicure_shared_modifiers& shmod);
