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

#pragma once

#include "starling_common/alignment.hh"
#include "starling_common/starling_read_segment.hh"
#include "starling_common/starling_shared.hh"

#include <vector>


/// for a read segment, expand into a set of starting alignments (exemplars)
/// for alignment search based on read mapper and grouper proposals
///
void
get_exemplar_alignments(const starling_options& opt,
                        const read_segment& rseg,
                        std::vector<alignment>& exal); 

