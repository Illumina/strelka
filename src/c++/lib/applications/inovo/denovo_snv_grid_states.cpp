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
