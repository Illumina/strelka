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

#ifndef __STARLING_TYPES_HH
#define __STARLING_TYPES_HH

#include <stdint.h>


typedef uint32_t align_id_t;
typedef int32_t sample_id_t;

// So long as the sample count required by multi-sample applications
// remains small, any typoe of heap allocation is a waste, so we use
// this to create stack arrays:
//
enum { MAX_SAMPLE=2 };

#endif
