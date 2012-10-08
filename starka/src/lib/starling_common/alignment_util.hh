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

#ifndef __ALIGNMENT_UTIL_HH
#define __ALIGNMENT_UTIL_HH


#include "blt_util/pos_range.hh"
#include "starling_common/alignment.hh"
#include "starling_common/indel.hh"


/// Provides the bounds of an alignment
///
/// note that for normal indel intersection test leading/trailing
/// indels will be missed if using this range
///
known_pos_range
get_strict_alignment_range(const alignment& al);

/// Provides the bounds of an alignment when edge indels
/// are converted to match.
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
/// start+read_length:
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

#if 0
/// Shift all indels as far "to the left" as possible -- note that
/// some indels may be lost. Returns true if the alignment was changed.
///
bool
normalize_alignment(alignment& al,
                    const std::string& read_seq,
                    const std::string& ref_seq);
#endif

#endif
