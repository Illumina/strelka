// -*- mode: c++; indent-tabs-mode: nil; -*-
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


#include "starling_common/IndelBuffer.hh"
#include "starling_common/starling_read.hh"
#include "starling_common/starling_base_shared.hh"
#include "ActiveRegionDetector.hh"


/// search for most likely realignments of the read and score alternate
/// indel states in preparation for indel genotype calling
///
/// Note search is for most likely realignment to a pool of potential
/// haplotypes formed by all possible candidate indel combinations. This
/// means there are not penalties for candidate indels if we think about the
/// realignment wrt the reference haplotype.
///
/// \param realign_buffer_range range in reference coordinates in which read is allowed to realign to (due to buffering constraints)
///
void
realign_and_score_read(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const starling_sample_options& sample_opt,
    const reference_contig_segment& ref,
    const known_pos_range& realign_buffer_range,
    const unsigned sampleId,
    const ActiveRegionDetector& activeRegionDetector,
    read_segment& rseg,
    IndelBuffer& indelBuffer);
