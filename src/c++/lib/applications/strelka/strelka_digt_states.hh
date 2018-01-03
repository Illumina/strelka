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

/// \file
/// \author Chris Saunders
/// \author Sangtae Kim
///

#pragma once

#include "blt_util/blt_types.hh"

#include <iosfwd>
#include <vector>


/// Enumerate allowed germline genotype states in the somatic model
namespace SOMATIC_DIGT
{
enum index_t
{
    REF,
    HOM,
    HET,
    SIZE
};

inline
const char*
label(const unsigned idx)
{
    switch (idx)
    {
    case REF:
        return "ref";
    case HOM:
        return "hom";
    case HET:
        return "het";
    default:
        return "xxx";
    }
}

}

namespace SOMATIC_STATE
{
enum index_t
{
    NON_SOMATIC,
    SOMATIC,
    SIZE
};

}

namespace DIGT_GRID
{

// HET_RES is the number of points sampled between 0 and 0.5 on
// the continuous frequency scale. Thus a fully sampled axis will
// be sampled HET_RES*2+3 times.
//

//// The first three states are designed to overlap with
//// STAR_DIINDEL (ie. conventional diploid indel model), after this
//// the grid frequency (ie. approximations of continuous frequency)
//// states are added. The grid states are treated just like the
//// STAR_DIINDEL het state for certain purposes (printing, for
//// instance)

enum constants { HET_RES = 9,
                 HET_COUNT = HET_RES*2+1,
                 HOM_SIZE = 2,
                 PRESTRAND_SIZE = HOM_SIZE+HET_COUNT,
                 STRAND_STATE_SIZE = HET_RES
               };

enum index_t { SIZE = PRESTRAND_SIZE+STRAND_STATE_SIZE };

const blt_float_t RATIO_INCREMENT = 0.5f/static_cast<blt_float_t>(DIGT_GRID::HET_RES+1);

blt_float_t
get_fraction_from_index(int index);

}

namespace DDIGT
{

enum index_t { SIZE = SOMATIC_DIGT::SIZE*SOMATIC_STATE::SIZE };

inline
unsigned
get_state(
    const unsigned normal_gt,
    const unsigned tumor_gt)
{
    return normal_gt*SOMATIC_STATE::SIZE + tumor_gt;
}

void
write_indel_state(const DDIGT::index_t dgt,
                  std::ostream& os);

void
write_snv_state(const DDIGT::index_t dgt,
                const char ref_base,
                const char normal_alt_base,
                const char tumor_alt_base,
                std::ostream& os);

inline
void
get_digt_states(
    const unsigned dgt,
    unsigned& normal_gt,
    unsigned& tumor_gt)
{
    normal_gt = (dgt/SOMATIC_STATE::SIZE);
    tumor_gt = (dgt%SOMATIC_STATE::SIZE);
}

void
write_alt_alleles(char normal_alt_base,
                  char tumor_alt_base,
                  char ref_base,
                  std::ostream& os);
}

namespace DDIGT_GRID
{

enum constants { PRESTRAND_SIZE = DIGT_GRID::PRESTRAND_SIZE*DIGT_GRID::PRESTRAND_SIZE };
enum index_t { SIZE = PRESTRAND_SIZE+DIGT_GRID::STRAND_STATE_SIZE };

inline
unsigned
get_state(
    const unsigned normal_freq,
    const unsigned tumor_freq)
{
    if (normal_freq<DIGT_GRID::PRESTRAND_SIZE) return normal_freq+DIGT_GRID::PRESTRAND_SIZE*tumor_freq;
    return PRESTRAND_SIZE+normal_freq-DIGT_GRID::PRESTRAND_SIZE;
}

struct is_nonsom_maker_t
{
    is_nonsom_maker_t();

    std::vector<bool> val;
};

extern const is_nonsom_maker_t is_nonsom;
}
