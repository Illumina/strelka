// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

/// \file
///
/// \author Chris Saunders
///

#ifndef __STARLING_READ_ALIGN_SCORE_HH
#define __STARLING_READ_ALIGN_SCORE_HH

#include "candidate_alignment.hh"

#include "starling_common/indel_buffer.hh"
#include "starling_common/starling_read.hh"
#include "starling_common/starling_shared.hh"



double
score_candidate_alignment(const starling_options& client_opt,
                          const indel_buffer& ibuff,
                          const read_segment& rseg,
                          const candidate_alignment& cal,
                          const reference_contig_segment& ref);

#endif
