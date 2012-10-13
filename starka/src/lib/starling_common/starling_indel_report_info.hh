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

#ifndef __STARLING_INDEL_REPORT_INFO_HH
#define __STARLING_INDEL_REPORT_INFO_HH

#include "starling_common/indel.hh"
#include "starling_common/starling_shared.hh"

#include <string>


// indel summary information which is shared between all samples:
//
struct starling_indel_report_info {

    starling_indel_report_info()
        : is_repeat_unit(false)
        , ref_repeat_count(0)
        , indel_repeat_count(0)
        , ihpol(0)
        , it(INDEL::NONE) {}

    std::string desc;
    std::string ref_seq;
    std::string indel_seq;
    std::string vcf_ref_seq;
    std::string vcf_indel_seq;
    std::string ref_upstream;
    std::string ref_downstream;

    bool is_repeat_unit;
    std::string repeat_unit;
    unsigned ref_repeat_count;
    unsigned indel_repeat_count;
    unsigned ihpol; // interupted homopolymer length

    // not directly reported, but handy to have pre-calculated:
    INDEL::index_t it;
};


// indel summary information which is specific to each sample:
//
struct starling_indel_sample_report_info {

    starling_indel_sample_report_info()
        : n_q30_ref_reads(0)
        , n_q30_indel_reads(0)
        , n_q30_alt_reads(0)
        , n_other_reads(0)
        , depth(0) {}

    unsigned n_q30_ref_reads;
    unsigned n_q30_indel_reads;
    unsigned n_q30_alt_reads;
    unsigned n_other_reads;
    unsigned depth;
};



// translate read path likelihoods to posterior probs
//
read_path_scores
indel_lnp_to_pprob(const starling_deriv_options& client_dopt,
                   const read_path_scores& read_path_lnp,
                   const bool is_tier2_pass,
                   const bool is_use_alt_indel);


/// get strings for desc, indel and ref sequences that span only the
/// indel (or indel combination) and other information used to
/// summarize indel output
///
void
get_starling_indel_report_info(const indel_key& ik,
                               const indel_data& id,
                               const reference_contig_segment& ref,
                               starling_indel_report_info& iri);


struct pos_basecall_buffer;

void
get_starling_indel_sample_report_info(const starling_deriv_options& dopt,
                                      const indel_key& ik,
                                      const indel_data& id,
                                      const pos_basecall_buffer& bc_buff,
                                      const bool is_include_tier2,
                                      const bool is_use_alt_indel,
                                      starling_indel_sample_report_info& isri);

#endif
