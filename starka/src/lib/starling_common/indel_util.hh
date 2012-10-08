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
