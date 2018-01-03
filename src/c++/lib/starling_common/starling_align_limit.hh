//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
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

#include <vector>


/// helps to control realignment complexity.
///
/// complexity is restricted as a function of the number of
/// indel candidates intersecting a read
///
struct starling_align_limit
{
    starling_align_limit(const unsigned max_alignment_count);

    unsigned
    get_max_toggle(const unsigned n_usable_indels) const
    {
        if (n_usable_indels>=_max_toggle.size())
        {
            return 1;
        }
        else
        {
            return _max_toggle[n_usable_indels];
        }
    }

private:
    std::vector<unsigned> _max_toggle;
};

