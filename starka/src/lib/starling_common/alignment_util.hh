// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

/// \file
///
/// \author Chris Saunders
///

#pragma once


#include "blt_util/pos_range.hh"
#include "starling_common/alignment.hh"
#include "starling_common/indel.hh"


/// Provides the bounds of an alignment
///
/// note that for normal indel intersection test leading/trailing
/// indels will be missed if using this range
///
/// (should this say insertions instead of indels?)
///
known_pos_range
get_strict_alignment_range(const alignment& al);

/// Provides the bounds of an alignment when edge indels
/// are converted to match, but not soft-clip segments
///
known_pos_range
get_soft_clip_alignment_range(const alignment& al);

/// Provides the bounds of an alignment when edge indels and soft clip
/// are converted to match.
///
known_pos_range
get_alignment_range(const alignment& al);


/// Provides the largest reasonable bounds of an alignment by
/// including any leading and trailing edge sequence and requiring that
/// the end of the alignment range is equal to at least the
/// start+read_length and the start of the alignment range is
/// equal to at least orig_end-read_length:
///
known_pos_range
get_alignment_zone(const alignment& al,
                   const unsigned seq_length);

/// Return true if the indel is in this alignment or is consistent
/// with an edge-insertion in the alignment.
///
/// read_indel_pr returns the left and right breakpoints in read
/// coordinates. The left or right value will not be set if the
/// breakpoint does not occur within the read.
///
bool
is_indel_in_alignment(const alignment& al,
                      const indel_key& ik,
                      pos_range& read_indel_pr);



/// remove any edge deletions and properly adjust pos for leading deletions:
alignment
remove_edge_deletions(const alignment& al,
                      const bool is_remove_leading_edge=true,
                      const bool is_remove_trailing_edge=true);


#if 0
/// Shift all indels as far "to the left" as possible -- note that
/// some indels may be lost. Returns true if the alignment was changed.
///
bool
normalize_alignment(alignment& al,
                    const std::string& read_seq,
                    const std::string& ref_seq);
#endif

