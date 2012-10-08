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

#ifndef __DEPTH_BUFFER_UTIL_HH
#define __DEPTH_BUFFER_UTIL_HH

#include "starling_common/alignment.hh"
#include "starling_common/depth_buffer.hh"


void
add_alignment_to_depth_buffer(const alignment& al,
                              depth_buffer& db);


#endif
