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

/// \file
/// \author Chris Saunders
///

#pragma once

#include "starling_common/indel_align_type.hh"
#include "starling_common/starling_pos_processor_base.hh"



/// Extract indel information from various alignment types and store
/// this information in the starling_pos_processor indel buffer.
///
/// \param edge_pin are the beginning or end of this read segment pinned?
///
/// assumes that path is already validated for read_seq!!!
///
unsigned
addAlignmentIndelsToPosProcessor(
    const unsigned max_indel_size,
    const reference_contig_segment& ref,
    const alignment& al,
    const bam_seq_base& read_seq,
    starling_pos_processor_base& posProcessor,
    const INDEL_ALIGN_TYPE::index_t iat,
    const align_id_t id,
    const unsigned sample_no,
    const std::pair<bool, bool>& edge_pin,
    const bool isLowMapQuality);

