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
/// \author Chris Saunders, Sangtae Kim
///

#include "strelka_digt_states.hh"

#include <iostream>



namespace DIGT_GRID
{

blt_float_t
get_fraction_from_index(int index)
{
    if (index == SOMATIC_DIGT::REF) return 0.f;
    if (index == SOMATIC_DIGT::HOM) return 1.f;
    if (index == SOMATIC_DIGT::HET) return 0.5f;
    if (index < SOMATIC_DIGT::SIZE+DIGT_GRID::HET_RES) return DIGT_GRID::RATIO_INCREMENT*(index-SOMATIC_DIGT::SIZE+1);
    return DIGT_GRID::RATIO_INCREMENT*(index-SOMATIC_DIGT::SIZE+2);
}

}

namespace DDIGT
{

void
write_indel_state(const DDIGT::index_t dgt,
                  std::ostream& os)
{
    unsigned normal_gt;
    unsigned tumor_gt;
    get_digt_states(dgt,normal_gt,tumor_gt);

    os << SOMATIC_DIGT::label(normal_gt);
    os << "->";
    if (tumor_gt == SOMATIC_STATE::NON_SOMATIC)
    {
        os << SOMATIC_DIGT::label(normal_gt);
    }
    else
    {
        if (normal_gt == SOMATIC_DIGT::REF)
            os << SOMATIC_DIGT::label(SOMATIC_DIGT::HET);   // ref->som is written as ref->het for backward compatibility
        else
            os << SOMATIC_DIGT::label(SOMATIC_DIGT::REF);   // het/hom->som is written as het/hom->ref for backward compatibility
    }
}

static
void
write_diploid_genotype(
    const char base1,
    const char base2,
    std::ostream& os)
{
    char diploid_genotype[3];
    if (base1 < base2)
    {
        diploid_genotype[0] = base1;
        diploid_genotype[1] = base2;
    }
    else
    {
        diploid_genotype[1] = base1;
        diploid_genotype[0] = base2;
    }
    diploid_genotype[2] = 0;    // null
    os << diploid_genotype;
}

void
write_snv_state(const DDIGT::index_t dgt,
                const char ref_base,
                const char normal_alt_base,
                const char tumor_alt_base,
                std::ostream& os)
{
    unsigned normal_gt;
    unsigned tumor_gt;
    get_digt_states(dgt,normal_gt,tumor_gt);

    switch (normal_gt)
    {
    case SOMATIC_DIGT::REF:
        write_diploid_genotype(ref_base, ref_base, os);
        os << "->";
        if (tumor_gt == SOMATIC_STATE::NON_SOMATIC)
            write_diploid_genotype(ref_base, ref_base, os);
        else
            write_diploid_genotype(ref_base, tumor_alt_base, os);
        break;
    case SOMATIC_DIGT::HET:
        write_diploid_genotype(ref_base, normal_alt_base, os);
        os << "->";
        if (tumor_gt == SOMATIC_STATE::NON_SOMATIC)
            write_diploid_genotype(ref_base, normal_alt_base, os);
        else
            write_diploid_genotype(ref_base, ref_base, os);
        break;
    case SOMATIC_DIGT::HOM:
        write_diploid_genotype(normal_alt_base, normal_alt_base, os);
        os << "->";
        if (tumor_gt == SOMATIC_STATE::NON_SOMATIC)
            write_diploid_genotype(normal_alt_base, normal_alt_base, os);
        else
            write_diploid_genotype(ref_base, ref_base, os);
        break;
    }
}

void
write_alt_alleles(char normal_alt_base,
                  char tumor_alt_base,
                  char ref_base,
                  std::ostream& os)
{
    if (tumor_alt_base != ref_base)
        // tumor: at least one non-ref call
        os << tumor_alt_base;
    else if (normal_alt_base != ref_base)
        // tumor: all reference calls, normal: at least one non-ref call
        os << normal_alt_base;
    else
        // tumor: all ref calls, normal: all ref calls
        os << ".";
}

}

namespace DDIGT_GRID
{
is_nonsom_maker_t::
is_nonsom_maker_t()
    : val(SIZE,false)
{
    for (unsigned gt(0); gt<DIGT_GRID::SIZE; ++gt)
    {
        val[get_state(gt,gt)] = true;
    }
}

const is_nonsom_maker_t is_nonsom;
}

