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

///
/// \author Chris Saunders
///

#pragma once

#include "blt_util/window_util.hh"


/////////////////////////////////
struct win_avg_set
{
    win_avg_set(const unsigned size)
        : ss_used_win(size)
        , ss_filt_win(size)
        , ss_spandel_win(size)
        , ss_submap_win(size)
    {}

    window_average ss_used_win;
    window_average ss_filt_win;
    window_average ss_spandel_win;
    window_average ss_submap_win;
};
