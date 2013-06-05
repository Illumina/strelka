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
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the beginning of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///

#pragma once


namespace INDEL_ALIGN_TYPE {
enum index_t {
    GENOME_TIER1_READ,
    GENOME_TIER2_READ,
    GENOME_SUBMAP_READ,
    CONTIG,
    CONTIG_READ
};

inline
const char*
label(const index_t i) {
    switch(i){
    case GENOME_TIER1_READ : return "genome_tier1";
    case GENOME_TIER2_READ : return "genome_tier2";
    case GENOME_SUBMAP_READ : return "genome_submap";
    case CONTIG : return "contig";
    case CONTIG_READ : return "contig_read";
    default: return "unknown";
    }
}
}

