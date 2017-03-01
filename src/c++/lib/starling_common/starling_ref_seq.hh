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

#include "starling_common/starling_base_shared.hh"
#include "blt_util/known_pos_range2.hh"
#include "htsapi/bam_header_info.hh"

#include <string>


void
setRefSegment(
    const starling_base_options& opt,
    const std::string& chrom,
    const known_pos_range2& range,
    reference_contig_segment& ref);


struct AnalysisRegionInfo
{
    /// chrom string from region parse
    std::string regionChrom;

    ///strelka range from region parse
    known_pos_range2 regionRange;

    ///analysis range padded by indel size (used for streamer classses)
    known_pos_range2 streamerRegionRange;

    ///analysis range padded by indel size + extra constant pad (used for reference region)
    known_pos_range2 refRegionRange;

    /// region formatted into a samtools region string, e.g. "chr20:100-200"
    std::string streamerRegion;
};


/// convert to string format expected by samtools/htslib
std::string
getSamtoolsRegionString(
    const std::string& chromName,
    const known_pos_range2& range);


/// given a genome segment for analysis described by chrom, beginPos, endPos,
/// produce various related region objects used to manage edge-effects
///
/// \param[in] beginPos start position (zero-indexed, closed)
/// \param[in] endPos end position (zero-indexed, open)
void
getStrelkaAnalysisRegionInfo(
    const std::string& chrom,
    const int32_t beginPos,
    const int32_t endPos,
    const unsigned maxIndelSize,
    AnalysisRegionInfo& rinfo);


/// given an input samtools region string for analysis, produce various related
/// region objects used to manage edge-effects
///
/// \param region[in] samtools formatted analysis region string
void
getStrelkaAnalysisRegionInfo(
    const std::string& region,
    const unsigned maxIndelSize,
    AnalysisRegionInfo& rinfo);


/// parse and sanity check regions
///
/// TODO reorg this into a different module
void
getStrelkaAnalysisRegions(
    const starling_base_options& opt,
    const std::string& referenceAlignmentFilename,
    const bam_header_info& referenceHeaderInfo,
    std::vector<AnalysisRegionInfo>& regionInfoList);
