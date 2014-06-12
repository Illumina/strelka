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
