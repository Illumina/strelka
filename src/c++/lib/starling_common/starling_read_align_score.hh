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

#include "candidate_alignment.hh"

#include "starling_common/indel_buffer.hh"
#include "starling_common/starling_read.hh"
#include "starling_common/starling_base_shared.hh"



double
score_candidate_alignment(const starling_base_options& client_opt,
                          const indel_buffer& ibuff,
                          const read_segment& rseg,
                          const candidate_alignment& cal,
                          const reference_contig_segment& ref);

