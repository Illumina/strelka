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

#include "CandidateAlignment.hh"

#include "starling_common/IndelBuffer.hh"
#include "starling_common/starling_read.hh"
#include "starling_common/starling_base_shared.hh"
#include "CandidateSnvBuffer.hh"

/// \return Score of candidate alignment \p cal for read segment \p rseg
///
/// The score is essentially `P(read | haplotype)`, where read=rseg and haplotype=ref+candidate alignment. This
/// is 'essentially' instead of 'exactly' because Strelka's haplotype model is not fully carried forward to
/// this point, so we have a collection of candidate indels and SNV positions rather than strict haplotype (set)
/// that the read is scored against. It is also an approximation because we're only considering a single alignment
/// of the read to the target haplotype.
///
double
scoreCandidateAlignment(
    const starling_base_options& opt,
    const IndelBuffer& indelBuffer,
    const read_segment& readSegment,
    const CandidateAlignment& cal,
    const reference_contig_segment& ref);

