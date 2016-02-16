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

#include <cstdint>
#include <vector>
#include <array>


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
    std::vector< std::array<uint8_t,2> > SampleGts;
    std::vector< std::string > gtstring;
    std::vector< unsigned > gqx;
    std::vector< unsigned > gq;
};
