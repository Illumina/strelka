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

/// \author Chris Saunders
///
#ifndef __MAP_LEVEL_HH
#define __MAP_LEVEL_HH

#include <cassert>


namespace MAPLEVEL {
enum index_t {
    UNKNOWN,
    UNKNOWN_MAPPED,
    TIER1_MAPPED,
    TIER2_MAPPED,
    SUB_MAPPED,
    UNMAPPED
};

inline
const char*
get_label(const index_t i) {
    switch (i) {
    case UNKNOWN: return "unknown";
    case UNKNOWN_MAPPED: return "unknown-mapped";
    case TIER1_MAPPED: return "tier1-mapped";
    case TIER2_MAPPED: return "tier2-mapped";
    case SUB_MAPPED: return "sub-mapped";
    case UNMAPPED: return "unmapped";
    default:
        assert(0);
        return "none";
    }
}
}


#endif
