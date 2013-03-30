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
