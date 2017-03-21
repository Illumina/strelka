//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
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

    void
    reset()
    {
        ss_used_win.reset();
        ss_filt_win.reset();
        ss_spandel_win.reset();
        ss_submap_win.reset();
    }

    window_average ss_used_win;
    window_average ss_filt_win;
    window_average ss_spandel_win;
    window_average ss_submap_win;
};
