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

    tier_t snv_tier = 0;
    tier_t snv_from_ntype_tier = 0;
    unsigned ref_gt = 0;
    result_set rs = result_set();
    bool is_forced_output = false;
    SiteNoise sn;
};
