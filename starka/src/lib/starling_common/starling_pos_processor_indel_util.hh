// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
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

#ifndef __STARLING_POS_PROCESSOR_INDEL_UTIL_HH
#define __STARLING_POS_PROCESSOR_INDEL_UTIL_HH

#include "starling_common/indel_align_type.hh"
#include "starling_common/starling_pos_processor_base.hh"


// Extract indel information from various alignment types and store
// this information in the starling_pos_processor indel buffer.
//
// assumes that path is already validated for read_seq!!!
//
void
add_alignment_indels_to_sppr(const unsigned max_indel_size,
                             const reference_contig_segment& ref,
                             const alignment& al,
                             const bam_seq_base& read_seq,
                             starling_pos_processor_base& sppr,
                             const INDEL_ALIGN_TYPE::index_t iat,
                             align_id_t id,
                             const unsigned sample_no,
                             const indel_set_t* edge_indel_ptr = NULL);

#endif
