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


#include "starling_common/indel_synchronizer.hh"
#include "starling_common/starling_read.hh"
#include "starling_common/starling_shared.hh"


///
/// \param realign_pr range in reference coordinates in which read is allowed to realign to (due to buffering constraints)
///
void
realign_and_score_read(const starling_options& opt,
                       const starling_deriv_options& dopt,
                       const starling_sample_options& sample_opt,
                       const reference_contig_segment& ref,
                       const known_pos_range& realign_pr,
                       read_segment& rseg,
                       indel_synchronizer& isync);

