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

#ifndef __STARLING_READ_ALIGN_CLIPPER_HH
#define __STARLING_READ_ALIGN_CLIPPER_HH

#include "candidate_alignment.hh"


typedef std::vector<const candidate_alignment*> cal_pool_t;


void
get_clipped_alignment_from_cal_pool(const cal_pool_t& max_cal_pool,
                                    const unsigned best_cal_id,
                                    alignment& al);

#endif
