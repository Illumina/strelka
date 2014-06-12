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

/// \file

/// variation on the original strowman snv caller -- implements a
/// compile-time specified grid in allele frequency space and requires
/// similar frequency as definition of non-somatic.
///

/// \author Chris Saunders
///

#pragma once


struct somatic_snv_genotype_grid {

    somatic_snv_genotype_grid()
        : is_forced_output(false) {}

    typedef bool tier_t;

    struct result_set {
        result_set()
            : snv_qphred(0)
            , snv_from_ntype_qphred(0)
            , nonsomatic_qphred(0)
        {}

        unsigned ntype;
        unsigned max_gt;
        int snv_qphred;
        int snv_from_ntype_qphred;
        int nonsomatic_qphred;
    };

    bool
    is_snv() const {
        return (0 != rs.snv_qphred);
    }

    bool
    is_output() const {
        return (is_snv() || is_forced_output);
    }

    tier_t snv_tier;
    tier_t snv_from_ntype_tier;
    unsigned ref_gt;
    result_set rs;
    bool is_forced_output;
};
