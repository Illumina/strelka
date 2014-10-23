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

/// variation on the original strawman snv caller -- implements a
/// compile-time specified grid in allele frequency space and requires
/// similar frequency as definition of non-somatic.
///

/// \author Chris Saunders
///

#pragma once

#include "SiteNoise.hh"


struct somatic_snv_genotype_grid
{
    typedef bool tier_t;

    struct result_set
    {
        unsigned ntype;
        unsigned max_gt;
        int snv_qphred = 0;
        int snv_from_ntype_qphred = 0;
        int nonsomatic_qphred = 0;

        double strandBias = 0;
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

    tier_t snv_tier;
    tier_t snv_from_ntype_tier;
    unsigned ref_gt;
    result_set rs;
    bool is_forced_output = false;
    SiteNoise sn;
};
