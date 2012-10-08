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
/// XXX_begin_pos is zero-indexed position at the begining of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///


#ifndef __STARLING_POS_PROCESSOR_CONTIG_UTIL_HH
#define __STARLING_POS_PROCESSOR_CONTIG_UTIL_HH

#include "blt_util/export_stream_reader.hh"
#include "starling_common/grouper_contig.hh"
#include "starling_common/starling_pos_processor_base.hh"

/// return true if contig can be used in sppr
bool
test_contig_usability(const starling_options& opt,
                      const grouper_contig& ctg,
                      const starling_pos_processor_base& sppr,
                      const char* sample_label = NULL);

/// pulls all contig reads corresponding to the current contig,
/// attempts to realign the read to the genome and if this succeeds,
/// adds the realigned read to the starling read buffer
///
/// assumes exr is already pointing to contig_id or a line that
/// procedes it, will skip any contig ids in between the current
/// position and target contig id.
///
void
process_contig_reads(const grouper_contig& ctg,
                     const unsigned max_indel_size,
                     export_stream_reader& exr,
                     starling_pos_processor_base& sppr,
                     bam_record& tmp_key_br,
                     const unsigned sample_no = 0);


void
process_contig(const starling_options& client_opt,
               const reference_contig_segment&,// ref,
               const grouper_contig& ctg,
               starling_pos_processor_base& sppr,
               const unsigned sample_no = 0,
               const char* sample_label = NULL);

#endif

