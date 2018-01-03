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

#include "HtsMergeStreamer.hh"
#include "htsapi/bam_record.hh"
#include "htsapi/bam_streamer.hh"
#include "htsapi/vcf_record.hh"
#include "options/AlignmentFileOptions.hh"
#include "starling_common/starling_pos_processor_base.hh"
#include "starling_common/starling_ref_seq.hh"

#include <functional>


/// \brief Register a set of alignment files to the hts streamer and verify consistency conditions.
std::vector<std::reference_wrapper<const bam_hdr_t> >
registerAlignments(
    const std::vector<std::string>& alignmentFilename,
    const std::vector<unsigned>& registrationIndices,
    HtsMergeStreamer& streamData);


/// re-segment each target region into 0-many target sub-regions, divided on each large gap
/// in bed-region coverage
///
/// this function assumes that call Bed region file is defined
///
/// \param[in] regionChrom chrom of the input genome segment
/// \param[in] regionRange range of the input genome segment
/// \param[out] subRegionRanges ranges for subregions of the input that clump together close bed file regions
///
void
getSubRegionsFromBedTrack(
    const std::string& callRegionsBedFilename,
    const std::string& regionChrom,
    const known_pos_range2& regionRange,
    std::vector<known_pos_range2>& subRegionRanges);


/// Handles input read alignments -- reads are parsed, their indels
/// are extracted and the reads/indels are buffered to posProcessor
///
void
processInputReadAlignment(
    const starling_base_options& opt,
    const reference_contig_segment& ref,
    const bam_streamer& read_stream,
    const bam_record& read,
    const pos_t base_pos,
    starling_read_counts& readCounts,
    starling_pos_processor_base& posProcessor,
    const unsigned sampleIndex = 0);


/// convert specified alt from vcf_record into an indel_observation
///
/// returns true if successful, false for exceeding max indel size or otherwise can't parse this VCF record
///
bool
convert_vcfrecord_to_indel_allele(
    const unsigned max_indel_size,
    const vcf_record& vcf_indel,
    const unsigned altIndex,
    IndelObservation& obs);


/// \brief Insert a candidate indel into posProcessor
///
/// This is typically used for candidate indel input from a vcf record.
///
/// \param[in] is_forced_output If true, then genotyping results must appear in the output for this indel.
void
process_candidate_indel(
    const unsigned max_indel_size,
    const vcf_record& vcf_indel,
    starling_pos_processor_base& posProcessor,
    const unsigned sampleIndex = 0,
    const bool is_forced_output = false);
