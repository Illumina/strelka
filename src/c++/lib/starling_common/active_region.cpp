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
/// \author Sangtae Kim
///

#include "active_region.hh"
#include <string>
#include <cstdio>

void active_region::insert_haplotype_base (align_id_t align_id, pos_t pos, std::string base)
{
    if (!_align_id_to_haplotype.count(align_id))
    {
        _align_id_to_haplotype[align_id] = std::string();
        for (int i=_start; i<pos; ++i)
            _align_id_to_haplotype[align_id] += '.';
    }
    _align_id_to_haplotype[align_id] += base;
    if (pos == _end)
        _align_id_reaching_end.insert(align_id);
}

void active_region::print_haplotypes() const
{
    std::map<std::string, unsigned> haplotype_counter;
    unsigned max_count = 0;
    for (auto entry : _align_id_to_haplotype)
    {
        align_id_t align_id = entry.first;
        std::string haplotype = entry.second;
        if (_align_id_reaching_end.find(align_id) == _align_id_reaching_end.end())
            haplotype += "*";

        if(!haplotype_counter.count(haplotype))
            haplotype_counter[haplotype] = 1;
        else
        {
            unsigned count = haplotype_counter[haplotype] + 1;
            if (count > max_count)
                max_count = count;
            haplotype_counter[haplotype] = count;
        }
    }

    for (auto entry : haplotype_counter)
    {
        std::string haplotype = entry.first;
        unsigned count = entry.second;
        if (count >= 3 && count >= max_count/4)
        {
            printf("%s\t%d\n", haplotype.c_str(), count);
        }
    }
}

