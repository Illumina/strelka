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


struct denovo_snv_call
{
    struct result_set
    {
        bool isCall = false;
    };

    bool
    is_snv() const
    {
        return rs.isCall;
    }

    bool
    is_output() const
    {
        return (is_snv() || is_forced_output);
    }

    unsigned ref_gt;
    result_set rs;
    bool is_forced_output = false;
};
