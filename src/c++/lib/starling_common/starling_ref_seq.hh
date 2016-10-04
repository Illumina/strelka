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

#include <string>


void
get_starling_ref_seq(
    const starling_base_options& opt,
    reference_contig_segment& ref);

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
    std::string streamerRegion;
};


/// given an input region for analysis, produce various related
/// region objects used to manage edge-effects
///
/// \param region[in] samtools formated analysis region string
void
getStrelkaAnalysisRegions(
    const std::string& region,
    const unsigned maxIndelSize,
    AnalysisRegionInfo& rinfo);
