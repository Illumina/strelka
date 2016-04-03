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

///
/// \author Chris Saunders
///

#pragma once

#include "alignment.hh"
#include "blt_util/reference_contig_segment.hh"
#include "htsapi/bam_record.hh"


/// Normalize alignment so that indels are left-shifted and reduced
///
/// Indels are left-shifted in such a way that the matching mismatch count is preserved. Multiple indels can
/// potentially be left-shifted into each other and joined
///
/// Indels are 'reduced' by identifying any combined insert/delete indels which can be matched to reference
///
/// Indels at or moved to the edge of an alignment are collapsed if possible, for indels on the left edge
/// of the alignment this could change the alignment position
///
/// All hard and soft clipping at the alignment edge is preserved.
///
/// Note that any BAM CIGAR seq-match/mismatch states ("=","X") will be collapsed to regular match ("M") states
/// in regions surrounding normalized indels
///
/// \returns true if the alignment is changed
bool
normalizeAlignment(
    const bam_seq_base& refSeq,
    const bam_seq_base& readSeq,
    alignment& al);

/// execute the above normalizeAlignment transformation directly on a bam_record
///
/// \returns true if the alignment is changed
bool
normalizeBamRecordAlignment(
    const reference_contig_segment& refSeq,
    bam_record& bamRead);
