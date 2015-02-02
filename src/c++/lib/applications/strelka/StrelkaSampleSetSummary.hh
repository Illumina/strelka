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

/// \file
///
/// \author Chris Saunders
///

///
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the beginning of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
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
