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

/// variation on the original strawman snv caller -- implements a
/// compile-time specified grid in allele frequency space and requires
/// similar frequency as definition of non-somatic.
///

/// \author Chris Saunders
///

#pragma once

#include "blt_util/digt.hh"
#include "blt_util/seq_util.hh"

#include <iosfwd>
#include <vector>


namespace DIGT_SGRID
{

// HET_RES is the number of points sampled between 0 and 0.5 on
// the continuous frequency scale. Thus a fully sampled axis will
// be sampled HET_RES*2+3 times.
//
// Note that the single-strand tests are only implemented on the
// half axis containing the reference allele as a major
// allele. The "STRAND_STATE_SIZE" appended to the PRESTAND_SIZE
// represents these additional strand-specific error states.
//
// Full layout:
//
// HET_STATE_SIZE - number of frequency levels:
// 10%, 20%, ...
//
// STRAND_STATE_SIZE - strand bias states:
// REF+1 10%, REF+2 10%, REF+3 10%
// REF+1 20%...
//


enum constants { HET_RES = 9,
                 HET_COUNT = HET_RES*2+1,
                 HOM_SIZE = 2,
                 PRESTRAND_SIZE = HOM_SIZE+HET_COUNT,
                 STRAND_STATE_SIZE = HET_RES
               };

enum index_t { SIZE = PRESTRAND_SIZE+STRAND_STATE_SIZE };

}

namespace DDIGT_SGRID
{

enum constants { PRESTRAND_SIZE = DIGT_SGRID::PRESTRAND_SIZE*DIGT_SGRID::PRESTRAND_SIZE };
enum index_t { SIZE = PRESTRAND_SIZE+DIGT_SGRID::STRAND_STATE_SIZE };

inline
unsigned
get_state(
    const unsigned normal_gt,
    const unsigned tumor_gt)
{
    if (normal_gt<DIGT_SGRID::PRESTRAND_SIZE) return normal_gt+DIGT_SGRID::PRESTRAND_SIZE*tumor_gt;
    return PRESTRAND_SIZE+normal_gt-DIGT_SGRID::PRESTRAND_SIZE;
}

inline
void
get_digt_grid_states(
    const unsigned dgt,
    unsigned& normal_gt,
    unsigned& tumor_gt)
{
    if (dgt<PRESTRAND_SIZE)
    {
        normal_gt = (dgt%DIGT_SGRID::PRESTRAND_SIZE);
        tumor_gt = (dgt/DIGT_SGRID::PRESTRAND_SIZE);
    }
    else
    {
        normal_gt= dgt+DIGT_SGRID::PRESTRAND_SIZE-PRESTRAND_SIZE;
        tumor_gt=normal_gt;
    }
}

void
write_state(const DDIGT_SGRID::index_t dgt,
            const unsigned ref_gt,
            std::ostream& os);

void
write_alt_alleles(unsigned alt_gt,
                  std::ostream& os);
}
