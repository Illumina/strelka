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

#include "blt_util/blt_types.hh"
#include "blt_util/reference_contig_segment.hh"


// Get the length of the longest homopolymer containing the current
// position if this position can be treated as any base.
//
unsigned
get_snp_hpol_size(const pos_t pos,
                  const reference_contig_segment& ref);

// find the largest homopolymer extending from pos where one
// occurance of an alternate base is allowed
//
unsigned
get_interupted_hpol_size(const pos_t pos,
                         const reference_contig_segment& ref);
