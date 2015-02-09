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


struct Tier2Options
{
    bool
    is_tier2() const
    {
        return
            (is_tier2_min_single_align_score ||
             is_tier2_min_paired_align_score ||
             is_tier2_single_align_score_rescue_mode ||
             is_tier2_mismatch_density_filter_count ||
             is_tier2_no_mismatch_density_filter ||
             is_tier2_no_filter_unanchored ||
             is_tier2_include_singleton ||
             is_tier2_include_anomalous ||
             is_tier2_indel_nonsite_match_prob);
    }

    int tier2_min_single_align_score = 0;
    bool is_tier2_min_single_align_score = false;
    int tier2_min_paired_align_score = 0;
    bool is_tier2_min_paired_align_score = false;
    bool is_tier2_single_align_score_rescue_mode = false;

    int tier2_mismatch_density_filter_count = 0;
    bool is_tier2_mismatch_density_filter_count = false;

    bool is_tier2_no_mismatch_density_filter = false;
    bool is_tier2_no_filter_unanchored = false;
    bool is_tier2_include_singleton = false;
    bool is_tier2_include_anomalous = false;

    bool is_tier2_indel_nonsite_match_prob = false;
    double tier2_indel_nonsite_match_prob = 0.25;
};
