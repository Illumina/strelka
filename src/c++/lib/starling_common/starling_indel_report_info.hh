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

#include "blt_common/MapqTracker.hh"
#include "starling_common/starling_base_shared.hh"
#include "starling_common/indel.hh"

#include <string>


namespace SimplifiedIndelReportType
{

/// some components of the indel reporting need to reduce all alternate alleles to the following
/// simplified states:
enum index_t
{
    INSERT,
    DELETE,
    SWAP,
    BREAKPOINT,
    OTHER
};

inline
index_t
getRateType(
    const IndelKey& indelKey)
{
    if     (indelKey.isPrimitiveDeletionAllele())
    {
        return DELETE;
    }
    else if (indelKey.isPrimitiveInsertionAllele())
    {
        return INSERT;
    }
    else if (indelKey.type == INDEL::INDEL)
    {
        return SWAP;
    }
    else if (indelKey.is_breakpoint())
    {
        return BREAKPOINT;
    }
    else
    {
        return OTHER;
    }
}
}


/// indel allele summary information which is shared between all samples:
///
struct starling_indel_report_info
{
    starling_indel_report_info() {}

    bool
    is_repeat_unit() const
    {
        return (! repeat_unit.empty());
    }

    void dump(std::ostream& os) const;

    std::string repeat_unit;
    unsigned repeat_unit_length = 0;
    unsigned ref_repeat_count = 0;
    unsigned indel_repeat_count = 0;
    unsigned ihpol = 0; ///< interrupted homopolymer length

    // not directly reported, but handy to have pre-calculated:
    SimplifiedIndelReportType::index_t it = SimplifiedIndelReportType::OTHER;
};

std::ostream& operator<<(std::ostream& os, const starling_indel_report_info& obj);


/// indel summary information which is specific to each sample:
///
struct starling_indel_sample_report_info
{
    /// TODO STREL-125 sample_report_info is still designed for only one alt allele

    starling_indel_sample_report_info() {}

    unsigned n_confident_ref_reads = 0;
    unsigned n_confident_indel_reads = 0;
    unsigned n_confident_alt_reads = 0;

    // number of lower-quality reads
    unsigned n_other_reads = 0;

    // the depth of the pileup preceding the indel
    unsigned tier1Depth = 0;

    // same as above, but by strand
    unsigned n_confident_ref_reads_fwd = 0;
    unsigned n_confident_indel_reads_fwd = 0;
    unsigned n_confident_alt_reads_fwd = 0;
    unsigned n_other_reads_fwd = 0;
    unsigned n_confident_ref_reads_rev = 0;
    unsigned n_confident_indel_reads_rev = 0;
    unsigned n_confident_alt_reads_rev = 0;
    unsigned n_other_reads_rev = 0;

    MapqTracker mapqTracker;

    ranksum readpos_ranksum;

    unsigned total_confident_reads() const
    {
        return n_confident_alt_reads + n_confident_indel_reads + n_confident_ref_reads;
    }

    void dump(std::ostream& os) const;
};
std::ostream& operator<<(std::ostream& os, const starling_indel_sample_report_info& obj);



/// translate read path likelihoods to posterior probs
///
ReadPathScores
indel_lnp_to_pprob(
    const starling_base_deriv_options& client_dopt,
    const ReadPathScores& read_path_lnp,
    const bool is_tier2_pass,
    const bool is_use_alt_indel);


/// get VCF specific ref and alt strings
void
getSingleIndelAlleleVcfSummaryStrings(
    const IndelKey& indelKey,
    const IndelData& indelData,
    const reference_contig_segment& ref,
    std::string& vcf_indel_seq,
    std::string& vcf_ref_seq);


/// get information used to summarize indel output
///
void
get_starling_indel_report_info(
    const IndelKey& indelKey,
    const reference_contig_segment& ref,
    starling_indel_report_info& indelReportInfo);


struct pos_basecall_buffer;

void
get_starling_indel_sample_report_info(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const IndelKey& indelKey,
    const IndelSampleData& indelSampleData,
    const pos_basecall_buffer& bc_buff,
    const bool is_include_tier2,
    const bool is_use_alt_indel,
    starling_indel_sample_report_info& isri);
