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

/// \file
///
/// \author Sangtae Kim
///

#include "active_region_detector.hh"
#include <cstdio>

void active_region_detector::insert_match(const align_id_t align_id, const pos_t pos, const char base_char)
{
    haplotype_base(align_id, pos) = std::string(1, base_char);
    add_align_id_to_pos(align_id, pos);
}

void
active_region_detector::insert_mismatch(const align_id_t align_id, const pos_t pos, const char base_char)
{
    add_count(pos, 1);
    haplotype_base(align_id, pos) = std::string(1, base_char);
    add_align_id_to_pos(align_id, pos);
}

void
active_region_detector::insert_indel(const indel_observation indel_obs)
{
    auto pos = indel_obs.key.pos;
    const int indel_count = 4;

    auto align_id = indel_obs.data.id;
    auto indel_key = indel_obs.key;
//    printf("indel\t%d\t%s\n", pos, indel_obs.data.insert_seq.c_str());
    if (indel_key.type == INDEL::INSERT)
    {
        add_count(pos-1, indel_count);
        add_count(pos, indel_count);
        haplotype_base(align_id, pos-1) += indel_obs.data.insert_seq;
//        printf("%d: insert %s to %d\n", align_id, indel_obs.data.insert_seq.c_str(), pos-1);
        add_align_id_to_pos(align_id, pos-1);
    }
    else if (indel_key.type == INDEL::DELETE)
    {
        unsigned length = indel_obs.key.length;
        for (unsigned i(0); i<length; ++i)
        {
            add_count(pos+i, indel_count);
            haplotype_base(align_id, pos+i) = "";
            add_align_id_to_pos(align_id, pos+i);
        }
        // length 1 deletion makes an active region
//        if (length == 1)
        add_count(pos-1, indel_count);
    }
    else
    {
        // ignore BP_LEFT, BP_RIGHT, SWAP
    }
}

//static int active_region_index = 0;

void
active_region_detector::update_start_position(const pos_t pos)
{
    for (int i(_buffer_start_pos); i<pos; ++i) reset_count(i);

    _buffer_start_pos = pos;

    if (_active_regions.empty()) return;

    active_region front =  _active_regions.front();
    if (front.get_start() == pos)
    {
//        printf("%d\t%d\t%d\t%d\t%d\n", ++active_region_index, front.get_start()+1, front.get_end()+1, front.get_num_variants(), front.get_length());
        front.print_haplotypes();
        _active_regions.pop_front();
    }
}

void
active_region_detector::update_end_position(const pos_t pos)
{
    bool is_current_pos_candidate_variant = is_candidate_variant(pos);
    if (pos - _active_region_start_pos  >= (int)_max_detection_window_size && pos - _prev_variant_pos > 1)
    {
        // this position doesn't extend the existing active region
        if (_num_variants >= _min_num_variants_per_region)
        {
            // close existing active region
            active_region activeRegion = active_region(_active_region_start_pos, _prev_variant_pos, _num_variants);
            _active_regions.push_back(activeRegion);
            // add haplotype bases
            printf(">%d\t%d\t%d\n", (_active_region_start_pos+1), (_prev_variant_pos+1), _num_variants);
            for (pos_t active_region_pos(_active_region_start_pos); active_region_pos<=_prev_variant_pos; ++active_region_pos)
            {
                for (align_id_t align_id : position_to_align_ids(active_region_pos))
                {
                    activeRegion.insert_haplotype_base(align_id, active_region_pos, haplotype_base(align_id, active_region_pos));
                }
            }
            activeRegion.print_haplotypes();
        }
        if (!is_current_pos_candidate_variant)
        {
            _active_region_start_pos = -1;
            _num_variants = 0;
        }
        else
        {
            // start new active region
            _active_region_start_pos = pos;
            _num_variants = 1;
        }
    }
    else
    {
        // this position doesn't extend the existing active region
        if (is_current_pos_candidate_variant)
        {
            ++_num_variants;
        }
    }

    if (is_current_pos_candidate_variant)
        _prev_variant_pos = pos;

    clear_pos(pos - _max_detection_window_size - 50);

    // max_deletion_size = 50
    // TODO: remaining
}

bool
active_region_detector::is_candidate_variant(const pos_t pos) const
{
    return get_count(pos) >= _min_num_mismatches_per_position;
}

std::list<active_region>&
active_region_detector::get_active_regions()
{
    return _active_regions;
}