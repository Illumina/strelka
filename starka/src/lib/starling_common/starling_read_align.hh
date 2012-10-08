// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
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
