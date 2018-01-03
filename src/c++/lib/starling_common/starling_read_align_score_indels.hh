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

///
/// \author Chris Saunders
///

#pragma once


#include "CandidateAlignment.hh"

#include "starling_common/IndelBuffer.hh"
#include "starling_common/starling_read_segment.hh"
#include "starling_common/starling_base_shared.hh"

#include <set>
#include <vector>


typedef std::map<IndelKey,bool> indel_status_map_t;


/// use the most likely alignment for each indel state for every indel
/// in indel_status_map to generate data needed in indel calling:
///
void
score_indels(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const starling_sample_options& sample_opt,
    const read_segment& rseg,
    IndelBuffer& indelBuffer,
    const unsigned sampleIndex,
    const std::set<CandidateAlignment>& candAlignments,
    const bool is_incomplete_search,
    const std::vector<double>& candAlignmentScores,
    double maxCandAlignmentScore,
    const CandidateAlignment* maxCandAlignmentPtr);
