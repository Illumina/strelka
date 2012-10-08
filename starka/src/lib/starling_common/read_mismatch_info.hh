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

#ifndef __READ_MISMATCH_INFO_HH
#define __READ_MISMATCH_INFO_HH


#include <vector>


struct read_base_mismatch_info {
    int mismatch_count;
    int mismatch_count_ns;

    // tmp workspace for create mismatch filter:
    int delta;

    bool mismatch_filter_map;
    bool tier2_mismatch_filter_map;

    // tmp workspace for create mismatch filter:
    bool is_mismatch;
};

typedef std::vector<read_base_mismatch_info> read_mismatch_info;


#endif
