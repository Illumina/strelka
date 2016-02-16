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

/// \author Chris Saunders
///

#include "position_somatic_snv_strand_grid_states.hh"
#include "blt_util/digt.hh"

#include <iostream>



namespace DIGT_SGRID
{

void
write_state(const DIGT_SGRID::index_t gt,
            const unsigned ref_gt,
            std::ostream& os)
{
    os << DIGT::label(DIGT_SGRID::get_digt_state(gt,ref_gt));
}

void
write_full_state(const DIGT_SGRID::index_t gt,
                 const unsigned ref_gt,
                 std::ostream& os)
{
    write_state(gt,ref_gt,os);
    os << "_" << get_het_count(gt);
    if (is_strand_state(gt)) os << "_strand";
}

strand_state_tables::
strand_state_tables()
{
    for (unsigned ref(0); ref<N_BASE; ++ref)
    {
        unsigned local_strand_state(0);
        for (unsigned alt(0); alt<N_BASE; ++alt)
        {
            if (alt==ref) continue;

            // search for the DIGT het state which contains
            // the given alt and ref allele:
            const unsigned gt(DIGT::get_gt_with_alleles(ref,alt));
            assert(gt>=N_BASE && gt<DIGT::SIZE);

            digt_state[ref][local_strand_state] = gt;
            strand_state[ref][gt-N_BASE] = local_strand_state;

            local_strand_state++;
        }
    }
}

const strand_state_tables stables;
}



namespace DDIGT_SGRID
{

void
write_state(const DDIGT_SGRID::index_t dgt,
            const unsigned /* ref_gt */,
            std::ostream& os)
{
    unsigned normal_gt;
    unsigned tumor_gt;
    DDIGT_SGRID::get_digt_grid_states(dgt,normal_gt,tumor_gt);

    os << DIGT_SIMPLE::label(normal_gt);
    os << "->";
    os << DIGT_SIMPLE::label(tumor_gt);
}

void
write_full_state(const DDIGT_SGRID::index_t dgt,
                 const unsigned ref_gt,
                 std::ostream& os)
{
    unsigned normal_gt;
    unsigned tumor_gt;
    DDIGT_SGRID::get_digt_grid_states(dgt,normal_gt,tumor_gt);

    DIGT_SGRID::write_full_state(static_cast<DIGT_SGRID::index_t>(normal_gt),ref_gt,os);
    os << "->";
    DIGT_SGRID::write_full_state(static_cast<DIGT_SGRID::index_t>(tumor_gt),ref_gt,os);
}


void
write_alt_alleles(unsigned alt_gt,
                  std::ostream& os)
{
    os << id_to_base(alt_gt);
}

is_nonsom_maker_t::
is_nonsom_maker_t()
    : val(SIZE,false)
{
    for (unsigned gt(0); gt<DIGT_SGRID::SIZE; ++gt)
    {
        val[get_state(gt,gt)] = true;
    }
}

const is_nonsom_maker_t is_nonsom;
}
