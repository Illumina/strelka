//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
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
/// \author Chris Saunders, Sangtae Kim
///


#pragma once

#include "boost/utility.hpp"
#include "SiteNoise.hh"

#include <iosfwd>

struct result_set
{
    unsigned ntype;
    unsigned max_gt;
    int qphred = 0;
    int from_ntype_qphred = 0;
};

std::ostream&
operator<<(std::ostream& os, const result_set& rs);

struct snv_result_set : result_set
{
    unsigned normal_alt_id;
    unsigned tumor_alt_id;
    int nonsomatic_qphred = 0;
    double strandBias = 0;
};

struct indel_result_set : result_set
{
    bool is_overlap = false;
};

struct somatic_snv_genotype_grid
{
    typedef bool tier_t;


    bool
    is_snv() const
    {
        return (0 != rs.qphred);
    }

    bool
    is_output() const
    {
        return (is_snv() || is_forced_output);
    }

    tier_t snv_tier = 0;
    tier_t snv_from_ntype_tier = 0;
    unsigned ref_gt = 0;
    snv_result_set rs = snv_result_set();
    bool is_forced_output = false;
    SiteNoise sn;
};

struct somatic_indel_call
{
    typedef bool tier_t;

    bool
    is_indel() const
    {
        return (rs.qphred != 0);
    }

    // should this indel be written out?
    bool
    is_output() const
    {
        return (is_indel() || is_forced_output);
    }

    tier_t sindel_tier;
    tier_t sindel_from_ntype_tier;
    indel_result_set rs;
    bool is_forced_output = false;
};
