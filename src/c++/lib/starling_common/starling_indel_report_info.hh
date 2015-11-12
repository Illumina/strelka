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

///
/// \author Chris Saunders
///

#pragma once

#include "starling_common/starling_base_shared.hh"
#include "starling_common/indel.hh"
#include <string>


/// indel summary information which is shared between all samples:
///
struct starling_indel_report_info
{
    std::string desc;
    std::string ref_seq;
    std::string indel_seq;
    std::string vcf_ref_seq;
    std::string vcf_indel_seq;
    std::string ref_upstream;
    std::string ref_downstream;

    bool
    is_repeat_unit() const
    {
        return (! repeat_unit.empty());
    }

    std::string repeat_unit;
    unsigned repeat_unit_length = 0;
    unsigned ref_repeat_count = 0;
    unsigned indel_repeat_count = 0;
    unsigned ihpol = 0; ///< interrupted homopolymer length

    // not directly reported, but handy to have pre-calculated:
    INDEL::index_t it = INDEL::NONE;

    void dump(std::ostream& os) const;

};
std::ostream& operator<<(std::ostream& os, const starling_indel_report_info& obj);


/// indel summary information which is specific to each sample:
///
struct starling_indel_sample_report_info
{
    unsigned n_q30_ref_reads = 0;
    unsigned n_q30_indel_reads = 0;
    unsigned n_q30_alt_reads = 0;

    // number of lower-quality reads
    unsigned n_other_reads = 0;

    // the depth of the pileup preceding the indel
    unsigned depth = 0;

    // same as above, but by strand
    unsigned n_q30_ref_reads_fwd = 0;
    unsigned n_q30_indel_reads_fwd = 0;
    unsigned n_q30_alt_reads_fwd = 0;
    unsigned n_other_reads_fwd = 0;
    unsigned n_q30_ref_reads_rev = 0;
    unsigned n_q30_indel_reads_rev = 0;
    unsigned n_q30_alt_reads_rev = 0;
    unsigned n_other_reads_rev = 0;

    double mean_mapq = 0.0;
    double mapq0_frac = 0.0;

    ranksum readpos_ranksum;

    unsigned total_q30_reads() const
    {
        return n_q30_alt_reads + n_q30_indel_reads + n_q30_ref_reads;
    }

    void dump(std::ostream& os) const;

};
std::ostream& operator<<(std::ostream& os, const starling_indel_sample_report_info& obj);



/// translate read path likelihoods to posterior probs
///
read_path_scores
indel_lnp_to_pprob(const starling_base_deriv_options& client_dopt,
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
get_starling_indel_sample_report_info(const starling_base_deriv_options& dopt,
                                      const indel_key& ik,
                                      const indel_data& id,
                                      const pos_basecall_buffer& bc_buff,
                                      const bool is_include_tier2,
                                      const bool is_use_alt_indel,
                                      starling_indel_sample_report_info& isri);
