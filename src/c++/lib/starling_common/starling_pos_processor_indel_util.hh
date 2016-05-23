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
/// \author Chris Saunders
///
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the beginning of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///

#pragma once

#include "starling_common/indel_align_type.hh"
#include "starling_common/starling_pos_processor_base.hh"


// Extract indel information from various alignment types and store
/// this information in the starling_pos_processor indel buffer.
///
/// \param edge_pin are the beginning or end of this read segment pinned?
///
/// assumes that path is already validated for read_seq!!!
///

struct read_stats
{

    read_stats(const uint8_t& map,
               const uint8_t* q)
//              const pos_t& p)
        : mapq(map),qual(q) {}

    const uint8_t& mapq;
    const uint8_t* qual;
//    const pos_t pos;
};



unsigned
add_alignment_indels_to_sppr(const unsigned max_indel_size,
                             const reference_contig_segment& ref,
                             const alignment& al,
                             const bam_seq_base& read_seq,
                             starling_pos_processor_base& sppr,
                             const INDEL_ALIGN_TYPE::index_t iat,
                             align_id_t id,
                             const unsigned sample_no,
                             const std::pair<bool,bool>& edge_pin,
                             const bool is_mapq_zero);

