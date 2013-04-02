// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

/// \file

/// \author Chris Saunders
///
#ifndef __SOMATIC_INDEL_CALL_HH
#define __SOMATIC_INDEL_CALL_HH

#include "boost/utility.hpp"


struct somatic_indel_call {

    somatic_indel_call()
        : is_indel(false)
    {}

    typedef bool tier_t;

    struct result_set {

        result_set()
            : is_overlap(false)
        {}

        unsigned ntype;
        unsigned max_gt;
        int sindel_qphred;
        int sindel_from_ntype_qphred;
        bool is_overlap;
    };

    bool is_indel;
    tier_t sindel_tier;
    tier_t sindel_from_ntype_tier;
    result_set rs;
};

#endif
