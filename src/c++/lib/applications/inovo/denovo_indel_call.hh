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

#include <cstdint>


namespace INOVO_TIERS
{
    enum index_t{
        TIER1,
        TIER2,
        SIZE
    };
}


struct denovo_indel_call
{
    struct result_set
    {
//        std::vector<unsigned> max_gt;
        unsigned max_gt; // this is the de-novo state, still need a mechanism to capture sample states
        int dindel_qphred = 0;
        bool is_overlap = false;
    };

    /// this is (one of) the criteria for writing out to a file rather than
    /// an indication the indel is PASS'd
    bool
    is_indel() const
    {
        return (rs.dindel_qphred != 0);
    }

    // should this indel be written out?
    bool
    is_output() const
    {
        return (is_indel() || is_forced_output);
    }

    result_set rs;
    uint8_t dindel_tier = 0;
    bool is_forced_output = false;
};
