// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
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
/// XXX_begin_pos is zero-indexed position at the begining of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///


#pragma once

#include "htsapi/bam_record.hh"
#include "htsapi/bam_streamer.hh"
#include "htsapi/vcf_record.hh"
#include "starling_common/starling_pos_processor_base.hh"


/// format the bam region string from the program options and
/// 'influence zones' defined by sppr.
///
std::string
get_starling_bam_region_string(const starling_options& opt,
                               const starling_deriv_options& dopt);


// handles genomic read alignments -- reads are parsed, their indels
// are extracted and buffered, and the reads themselves are buffered
//
void
process_genomic_read(const starling_options& client_opt,
                     const reference_contig_segment&, // ref,
                     const bam_streamer& read_stream,
                     const bam_record& read,
                     const pos_t base_pos,
                     const pos_t,// report_begin_pos,
                     starling_read_counts&, // brc,
                     starling_pos_processor_base& sppr,
                     const unsigned sample_no = 0);

/// insert a candidate indel into sppr
///
/// \param is_forced_output - the results of the genotype type must be output for this indel, no matter how unlikely the variant is:
void
process_candidate_indel(
    const unsigned max_indel_size,
    const vcf_record& vcf_indel,
    starling_pos_processor_base& sppr,
    const unsigned sample_no = 0,
    const bool is_forced_output = false);
