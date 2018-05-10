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

#include "starling_common/SampleSetSummary.hh"


namespace STRELKA_SAMPLE_TYPE
{
enum index_t { NORMAL, TUMOR, SIZE };

inline
const char*
get_label(const unsigned i)
{
    switch (static_cast<index_t>(i))
    {
    case NORMAL:
        return "NORMAL";
    case TUMOR:
        return "TUMOR";
    default:
        return "UNKNOWN";
    }
}

inline
char
get_char_label(const unsigned i)
{
    switch (static_cast<index_t>(i))
    {
    case NORMAL:
        return 'n';
    case TUMOR:
        return 't';
    default:
        return '?';
    }
}
}


// same thing, but easier to pass around as an argument:
//
struct StrelkaSampleSetSummary : public SampleSetSummary
{
    StrelkaSampleSetSummary() : SampleSetSummary() {}

    unsigned
    size() const override
    {
        return STRELKA_SAMPLE_TYPE::SIZE;
    }

    const char*
    get_label(
        const unsigned i) const override
    {
        return STRELKA_SAMPLE_TYPE::get_label(i);
    }

    const char*
    get_prefix(
        const unsigned i,
        const bool is_tier1) const override
    {
        using namespace STRELKA_SAMPLE_TYPE;

        switch (static_cast<index_t>(i))
        {
        case NORMAL:
            return (is_tier1 ? "n1-" : "n2-");
        case TUMOR:
            return (is_tier1 ? "t1-" : "t2-");
        default:
            return "?" "?-";
        }
    }
};
