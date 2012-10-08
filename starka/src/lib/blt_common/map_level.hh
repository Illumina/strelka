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
        switch(i) {
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
