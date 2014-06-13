// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
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

#ifndef __READ_MISMATCH_INFO_HH
#define __READ_MISMATCH_INFO_HH


#include <vector>


struct read_base_mismatch_info
{
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
