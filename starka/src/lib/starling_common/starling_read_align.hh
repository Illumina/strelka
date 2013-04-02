// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

/// \file
///
/// \author Chris Saunders
///

#ifndef __STARLING_READ_ALIGN_HH
#define __STARLING_READ_ALIGN_HH

#include "starling_common/indel_synchronizer.hh"
#include "starling_common/starling_read.hh"
#include "starling_common/starling_shared.hh"

#include <string>


void
realign_and_score_read(const starling_options& opt,
                       const starling_deriv_options& dopt,
                       const starling_sample_options& sample_opt,
                       const reference_contig_segment& ref,
                       read_segment& rseg,
                       indel_synchronizer& isync);

#endif
