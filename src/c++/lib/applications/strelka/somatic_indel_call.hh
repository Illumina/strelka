// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

/// \author Chris Saunders
///

#pragma once

#include "boost/utility.hpp"


struct somatic_indel_call {

    somatic_indel_call()
        : is_forced_output(false)
    {}

    typedef bool tier_t;

    struct result_set {

        result_set()
            : sindel_qphred(0),
              is_overlap(false)
        {}

        unsigned ntype;
        unsigned max_gt;
        int sindel_qphred;
        int sindel_from_ntype_qphred;
        bool is_overlap;
    };

    bool
    is_indel() const {
        return (rs.sindel_qphred != 0);
    }

    // should this indel be written out?
    bool
    is_output() const {
        return (is_indel() || is_forced_output);
    }

    tier_t sindel_tier;
    tier_t sindel_from_ntype_tier;
    result_set rs;
    bool is_forced_output;
};
