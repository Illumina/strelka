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

#include "htsapi/bam_record.hh"


namespace READ_FILTER_TYPE
{
enum index_t
{
    PRIMARY,
    DUPLICATE,
    UNMAPPED,
    SECONDARY,
    SUPPLEMENTARY,
    NONE
};

inline
const char*
label(const index_t id)
{
    switch (id)
    {
    case PRIMARY:
        return "Primary";
    case DUPLICATE:
        return "Duplicate";
    case UNMAPPED:
        return "Unmapped";
    case SECONDARY:
        return "Secondary";
    case SUPPLEMENTARY:
        return "Supplementary";
    default:
        return "None";
    }
}
}


/// read filters which are *always* used, because starling/strelka
/// can't do anything sensible with this information:
///
inline
READ_FILTER_TYPE::index_t
starling_read_filter_shared(
    const bam_record& bamRead)
{
    using namespace READ_FILTER_TYPE;

    if (bamRead.is_filter()) return PRIMARY;
    if (bamRead.is_dup()) return DUPLICATE;
    if (bamRead.is_unmapped()) return UNMAPPED;
    if (bamRead.is_secondary()) return SECONDARY;
    if (bamRead.is_supplementary()) return SUPPLEMENTARY;

    return NONE;
}
