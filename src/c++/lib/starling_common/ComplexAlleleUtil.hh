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
/// \author Sangtae Kim
///

#pragma once

#include "ActiveRegion.hh"
#include "starling_pos_processor_base.hh"

unsigned addComplexAlleleToSppr(const reference_contig_segment &ref,
                                const alignment &al,
                                const bam_seq_base &read_seq,
                                starling_pos_processor_base& sppr,
                                const INDEL_ALIGN_TYPE::index_t iat,
                                const align_id_t id,
                                const unsigned sample_no,
                                const std::pair<bool, bool> &edge_pin,
                                const bool is_mapq_zero);

