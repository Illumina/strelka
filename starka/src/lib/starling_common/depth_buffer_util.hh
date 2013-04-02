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

#ifndef __DEPTH_BUFFER_UTIL_HH
#define __DEPTH_BUFFER_UTIL_HH

#include "starling_common/alignment.hh"
#include "starling_common/depth_buffer.hh"


void
add_alignment_to_depth_buffer(const alignment& al,
                              depth_buffer& db);


#endif
