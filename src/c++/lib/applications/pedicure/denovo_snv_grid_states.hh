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

/// facilitates a compile-time specified grid in allele frequency space
///

/// \author Chris Saunders
///

#pragma once

#include "blt_util/digt.hh"
#include "blt_util/seq_util.hh"

#include <iosfwd>
#include <vector>


namespace DIGT_DGRID
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
// N_BASE - homozygous genotypes:
// AA, CC, GG, TT
//
// HET_STATE_SIZE - number of het types x number of frequency levels:
// AC 10%, AG 10%, AT 10%, CG 10%...
// AC 20%...
//
// STRAND_STATE_SIZE - strand bias states:
// REF+1 10%, REF+2 10%, REF+3 10%
// REF+1 20%...
//


enum constants { HET_RES = 4,
                 HET_COUNT = HET_RES*2+1,
                 HET_STATE_SIZE = DIGT::HET_SIZE*HET_COUNT,
               };

enum index_t { SIZE = N_BASE+HET_STATE_SIZE };



/// return mixture frequency index (ie. 0= 10% 1= 20%, etc...)
inline
unsigned
get_het_count(const unsigned state)
{
    assert(state<SIZE);
    if (state<N_BASE) return 0;
    return (state-N_BASE)/DIGT::HET_SIZE;
}

/// reduce complex strelka state to closest diploid state based on non-zero alleles present in the strelka state
inline
unsigned
get_digt_state(const unsigned state)
{
    assert(state<SIZE);
    if (state<N_BASE)         return state;
    return N_BASE+((state-N_BASE)%DIGT::HET_SIZE);
}

// write only most representative genotype, ie "GT", for each state
void
write_state(
    const DIGT_DGRID::index_t gt,
    std::ostream& os);

// write genotype, heterozygous id
void
write_full_state(
    const DIGT_DGRID::index_t gt,
    std::ostream& os);

}
