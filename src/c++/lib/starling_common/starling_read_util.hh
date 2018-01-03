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

///
/// \author Chris Saunders
///

#pragma once

#include "blt_common/blt_shared.hh"
#include "blt_util/blt_types.hh"
#include "htsapi/bam_seq.hh"
#include "starling_common/alignment.hh"
#include "starling_common/read_mismatch_info.hh"
#include "ActiveRegionDetector.hh"


// the position used to buffer alignments allows for the possibility that
// any leading soft-clip or insert sequence could really align to the reference:
//
pos_t
get_alignment_buffer_pos(const alignment& al);


/// \brief generate mismatch-filter map in a way that scales
/// efficiently for large mismatch density windows:
///
void
create_mismatch_filter_map(const blt_options& client_opt,
                           const alignment& al,
                           const bam_seq_base& ref_seq,
                           const bam_seq_base& read_seq,
                           const unsigned read_begin,
                           const unsigned read_end,
                           const CandidateSnvBuffer& candidateSnvBuffer,
                           read_mismatch_info& rmi);


/// \brief find the range within this alignment which is considered
/// sufficiently high quality to support indel candidacy.
///
/// Note this function only guards against random/noise alignments,
/// sub-optimal alignments with sufficient supporting sequence on each
/// side is not identified.
///
/// if nothing is excluded, then valid_pr.begin_pos=0 and
/// valid_pr.end_pos=read_size
///
void
get_valid_alignment_range(const alignment& al,
                          const bam_seq_base& ref_seq,
                          const bam_seq_base& read_seq,
                          pos_range& valid_pr);
