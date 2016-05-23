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

#pragma once

#include "blt_util/blt_types.hh"

class active_region
{
public:
    active_region(pos_t start, pos_t end, unsigned num_variants):
        _start(start), _end(end), _num_variants(num_variants)
    {}

    pos_t get_start() { return _start; }
    pos_t get_end() { return _end; }
    unsigned get_length() { return _end - _start + 1; }
    unsigned get_num_variants() { return _num_variants; }
    bool contains(pos_t pos) { return pos >= _start && pos <= _end; }

private:
    pos_t _start;
    pos_t _end;
    unsigned _num_variants;
};
