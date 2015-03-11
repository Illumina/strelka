// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
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

