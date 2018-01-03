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

#include "alignment.hh"
#include "blt_util/reference_contig_segment.hh"
#include "htsapi/bam_record.hh"


/// Normalize alignment so that indels are left-shifted and simplified.
///
/// Indels are left-shifted in such a way that the match count of the alignment
/// does not get worse with each 1-base left-shift operation.
///
/// Multiple indels can potentially be left-shifted to become adjacent. If this
/// happens the indels will be merged into one larger indel.
///
/// Indels are 'simplified' by identifying any portion of a set of adjacent
/// insert/delete edits which can be matched to reference. For example:
/// REF: ACTGC, READ: ACGC, CIGAR: 2M1I2D1M
/// Should collapse to CIGAR: 2M1D2M
///
/// Indels at or left-shifted to the edge of an alignment are simplified if possible.
/// For indels on the left edge of the alignment this could change the alignment
/// position.
///
/// All hard and soft clipping at the alignment edge is preserved.
///
/// Note that any BAM CIGAR seq-match/mismatch states ("=","X") will be collapsed to
/// regular match ("M") states in regions surrounding normalized indels.
///
/// \returns True if the input alignment is changed
///
bool
normalizeAlignment(
    const bam_seq_base& refSeq,
    const bam_seq_base& readSeq,
    alignment& al);

/// Execute the above normalizeAlignment transformation directly on a bam_record
///
/// \returns True if bamRead's alignment is changed
///
bool
normalizeBamRecordAlignment(
    const reference_contig_segment& refSeq,
    bam_record& bamRead);
