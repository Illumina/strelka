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

///
/// \author Chris Saunders, Sangtae Kim
///

#include "strelka_digt_states.hh"

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
            os << SOMATIC_DIGT::label(SOMATIC_DIGT::HET);   // ref->som is written as ref->het for backward compatability
        else
            os << SOMATIC_DIGT::label(SOMATIC_DIGT::REF);   // het/hom->som is written as het/hom->ref for backward compatability
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

    switch (normal_gt) {
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
write_alt_alleles(unsigned alt_gt,
                  std::ostream& os)
{
    os << id_to_base(alt_gt);
}

}
