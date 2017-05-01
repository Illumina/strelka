//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
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

void
getAlleleSampleReportInfo(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const IndelKey& indelKey,
    const IndelSampleData& indelSampleData,
    const pos_basecall_buffer& bc_buff,
    const bool is_include_tier2,
    const bool is_use_alt_indel,
    AlleleSampleReportInfo& isri);


unsigned
shortestUnencountered(
		      const reference_contig_segment& ref,
		      const pos_t pos,
		      const unsigned numEncoded,
		      const bool left = false);

unsigned
computeContextCompressability(
			      const reference_contig_segment& ref,
			      const pos_t leftPos,
			      const pos_t rightPos,
			      const unsigned numKeys);
