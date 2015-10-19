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
/// \author Morten Kallberg
///

#pragma once

#include <cstdint>
#include <array>
#include <vector>


struct denovo_snv_call
{
    struct result_set
    {
        unsigned max_gt;
        int dsnv_qphred = 0;
        int snv_qphred = 1;
    };

    bool
    is_snv() const
    {
        return (0 != rs.snv_qphred);
    }

    bool
    is_output() const
    {
        return (is_snv() || is_forced_output);
    }

    unsigned ref_gt = 0;
    uint8_t dsnv_tier = 0;
    bool is_forced_output = false;
    result_set rs;
    
    std::vector< std::array<float,3> > Sampleplhoods;

};
