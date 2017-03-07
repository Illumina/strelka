// -*- mode: c++; indent-tabs-mode: nil; -*-
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

/// \author Chris Saunders
///

#include "denovo_snv_grid_states.hh"

#include "blt_util/digt.hh"

#include <iostream>



namespace DIGT_DGRID
{

void
write_state(const DIGT_DGRID::index_t gt,
            std::ostream& os)
{
    os << DIGT::label(DIGT_DGRID::get_digt_state(gt));
}



void
write_full_state(const DIGT_DGRID::index_t gt,
                 std::ostream& os)
{
    write_state(gt,os);
    os << "_" << get_het_count(gt);
}

}
