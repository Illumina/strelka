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
                             const indel_set_t* edge_indel_ptr = NULL);

