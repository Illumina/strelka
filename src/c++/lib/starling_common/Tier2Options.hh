//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
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

#pragma once


struct Tier2Options
{
    bool
    is_tier2() const
    {
        return
            (is_tier2_min_mapping_quality ||
             is_tier2_mismatch_density_filter_count ||
             is_tier2_no_mismatch_density_filter ||
             is_tier2_include_singleton ||
             is_tier2_include_anomalous ||
             is_tier2_indel_nonsite_match_prob);
    }

    int tier2_min_mapping_quality = 0;
    bool is_tier2_min_mapping_quality = false;

    int tier2_mismatch_density_filter_count = 0;
    bool is_tier2_mismatch_density_filter_count = false;

    bool is_tier2_no_mismatch_density_filter = false;
    bool is_tier2_include_singleton = false;
    bool is_tier2_include_anomalous = false;

    bool is_tier2_indel_nonsite_match_prob = false;
    double tier2_indel_nonsite_match_prob = 0.25;
};
