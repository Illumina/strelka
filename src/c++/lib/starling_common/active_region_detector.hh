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

#pragma once

#include "active_region.hh"
#include "blt_util/blt_types.hh"
#include "starling_read_segment.hh"
#include "indel.hh"
#include <vector>
#include <list>
#include <set>

class active_region_detector {
public:
    static const unsigned max_buffer_size = 1000;
    static const unsigned max_depth = 1000;

    active_region_detector(
            unsigned max_detection_window_size = 30,
            unsigned min_num_mismatches_per_position = 9,
            unsigned min_num_variants_per_region = 2) :
            _max_detection_window_size(max_detection_window_size),
            _min_num_mismatches_per_position(min_num_mismatches_per_position),
            _min_num_variants_per_region(min_num_variants_per_region),
            _active_regions(),
            _variant_counter(max_buffer_size),
            _align_ids_current_active_region(),
            _position_to_align_ids(max_buffer_size, std::set<align_id_t>()),
            _haplotype_base(max_depth, std::vector<std::string>(max_buffer_size, std::string()))

    {
        _buffer_start_pos = 0;

        _num_variants = 0;
        _active_region_start_pos = -max_buffer_size;
        _prev_variant_pos = -_max_detection_window_size - 1;
    }

    void insert_match(const align_id_t align_id, const pos_t pos, const char base_code);
    void insert_mismatch(const align_id_t align_id, const pos_t pos, const char base_char);
    void insert_indel(const indel_observation indel_obs);
    void update_start_position(const pos_t pos);
    void update_end_position(const pos_t pos);
    bool is_empty() const { return _active_regions.empty(); }
    std::list<active_region> &get_active_regions();

private:
    unsigned _max_detection_window_size;
    unsigned _min_num_mismatches_per_position;
    unsigned _min_num_variants_per_region;

    pos_t _buffer_start_pos;

    pos_t _prev_variant_pos;
    pos_t _active_region_start_pos;
    unsigned _num_variants;

    std::list<active_region> _active_regions;
    std::vector<unsigned> _variant_counter;

    // for haplotypes
    std::vector<align_id_t> _align_ids_current_active_region;
    std::vector<std::set<align_id_t>> _position_to_align_ids;
    std::vector<std::vector<std::string>> _haplotype_base;

    bool is_candidate_variant(const pos_t pos) const;

    inline void reset_count(const pos_t pos)
    {
        _variant_counter[pos % max_buffer_size] = 0;
    }

    inline void add_count(const pos_t pos, unsigned count = 1)
    {
        _variant_counter[pos % max_buffer_size] += count;
    }

    inline unsigned get_count(const pos_t pos) const
    {
        return _variant_counter[pos % max_buffer_size];
    }

    inline void add_align_id_to_pos(const align_id_t align_id, const pos_t pos)
    {
        _position_to_align_ids[pos % max_buffer_size].insert(align_id);
    }

    inline std::set<align_id_t>& position_to_align_ids(const pos_t pos)
    {
        return _position_to_align_ids[pos % max_buffer_size];
    }

    inline std::string& haplotype_base(const align_id_t id, const pos_t pos)
    {
        return _haplotype_base[id % max_depth][pos % max_buffer_size];
    }

    inline void clear_pos(pos_t pos)
    {
        _position_to_align_ids[pos % max_buffer_size].clear();
    }
};



