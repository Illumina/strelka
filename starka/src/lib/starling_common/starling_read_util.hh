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

#ifndef __STARLING_READ_UTIL_HH
#define __STARLING_READ_UTIL_HH


#include "blt_common/blt_shared.hh"
#include "blt_util/bam_seq.hh"
#include "blt_util/blt_types.hh"
#include "blt_util/export_line_parser.hh"
#include "starling_common/alignment.hh"
#include "starling_common/read_mismatch_info.hh"

#include <string>


// provide a key from the export line which follows the pattern used
// in export->SAM/BAM conversion
//
void
get_read_key_from_export_line(const export_line_parser& exl,
                              std::string& key);

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
                           read_mismatch_info& rmi);


/// \brief find the range within this alignment which is considered
/// sufficiently high quality to support indel candidacy.
///
/// Note this function only gaurds against random/noise alignments,
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

#endif
