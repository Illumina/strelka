// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
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

/// \file
///
/// \author Chris Saunders
///

#pragma once

#include "starling_common/alignment.hh"
#include "starling_common/indel_set.hh"


struct candidate_alignment
{
    bool
    operator<(const candidate_alignment& rhs) const
    {
        if (al < rhs.al) return true;
        if (al == rhs.al)
        {
            if (leading_indel_key < rhs.leading_indel_key) return true;
            if (leading_indel_key == rhs.leading_indel_key)
            {
                return (trailing_indel_key < rhs.trailing_indel_key);
            }
        }
        return false;
    }

    alignment al;
    indel_key leading_indel_key;
    indel_key trailing_indel_key;
};


std::ostream& operator<<(std::ostream& os, const candidate_alignment& cal);



// get the keys of the indels present in the candidate alignment
//
void
get_alignment_indels(const candidate_alignment& cal,
                     const unsigned max_indel_size,
                     indel_set_t& indels);
