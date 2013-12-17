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

#ifndef __STARLING_READ_ALIGN_CLIPPER_HH
#define __STARLING_READ_ALIGN_CLIPPER_HH

#include "candidate_alignment.hh"


typedef std::vector<const candidate_alignment*> cal_pool_t;


void
get_clipped_alignment_from_cal_pool(const cal_pool_t& max_cal_pool,
                                    const unsigned best_cal_id,
                                    alignment& al);

#endif
