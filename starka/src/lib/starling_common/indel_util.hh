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

#ifndef __INDEL_UTIL_HH
#define __INDEL_UTIL_HH

#include "blt_util/pos_range.hh"
#include "starling_common/indel.hh"


// do these two indels overlap?:
bool
is_indel_conflict(const indel_key& ik1,
                  const indel_key& ik2);

bool
is_range_intersect_indel_breakpoints(const known_pos_range read_pr,
                                     const indel_key& ik);

#endif
