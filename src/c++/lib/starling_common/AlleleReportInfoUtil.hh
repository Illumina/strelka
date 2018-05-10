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

#include "AlleleReportInfo.hh"
#include "indel.hh"
#include "blt_util/reference_contig_segment.hh"
#include "starling_common/starling_base_shared.hh"


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
getAlleleReportInfo(
    const IndelKey& indelKey,
    const reference_contig_segment& ref,
    AlleleReportInfo& indelReportInfo);


struct pos_basecall_buffer;

/// Compute the (now legacy) AlleleSampleReportInfo data structure values for a specific indel and sample
///
/// The principal reason this structure is outdated is that it is designed in an allele, rather than locus
/// specific fashion. It does not fully support overlapping alleles directly (although there is some indirect
/// support)
///
/// \param[in] isUseTier2Data If true, include tier2 data while updating indel posteriors, otherwise only tier1 is used
/// \param[in] isUseAltIndel If true, account for indels other than indelKey while computing indel posteriors
void
getAlleleSampleReportInfo(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const IndelKey& indelKey,
    const IndelSampleData& indelSampleData,
    const pos_basecall_buffer& bc_buff,
    const bool isUseTier2Data,
    const bool isUseAltIndel,
    AlleleSampleReportInfo& isri);
