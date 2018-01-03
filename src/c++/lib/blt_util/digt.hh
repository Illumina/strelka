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
/// \author Chris Saunders
///

#pragma once

#include "seq_util.hh"

#include <cassert>

namespace DIGT
{
enum index_t
{
    AA,
    CC,
    GG,
    TT,
    AC,
    AG,
    AT,
    CG,
    CT,
    GT,
    SIZE
};

enum constants
{
    HET_SIZE = SIZE-N_BASE,
};

inline
const char*
label(const unsigned idx)
{
    switch (idx)
    {
    case AA:
        return "AA";
    case CC:
        return "CC";
    case GG:
        return "GG";
    case TT:
        return "TT";
    case AC:
        return "AC";
    case AG:
        return "AG";
    case AT:
        return "AT";
    case CG:
        return "CG";
    case CT:
        return "CT";
    case GT:
        return "GT";
    default:
        return "XX";
    }
}

inline
bool
is_het(const unsigned idx)
{
    return (idx>=N_BASE);
}

inline
char
hom_base(const unsigned idx)
{
    if (is_het(idx)) return 'X';
    return id_to_base(idx);
}

/// provide the alelle frequency for base_id given a diploid genotype gt
inline
double
expect(
    const int base_id,
    const int gt)
{
    static const double ex[SIZE][N_BASE] =
    {
        { 1.0, 0.0, 0.0, 0.0},
        { 0.0, 1.0, 0.0, 0.0},
        { 0.0, 0.0, 1.0, 0.0},
        { 0.0, 0.0, 0.0, 1.0},
        { 0.5, 0.5, 0.0, 0.0},
        { 0.5, 0.0, 0.5, 0.0},
        { 0.5, 0.0, 0.0, 0.5},
        { 0.0, 0.5, 0.5, 0.0},
        { 0.0, 0.5, 0.0, 0.5},
        { 0.0, 0.0, 0.5, 0.5}
    };

    return ex[gt][base_id];
}

/// coded form of expect function above:
/// 0 is 0.0
/// 1 is 0.5
/// 2 is 1.0
inline
unsigned
expect2(const int base_id,
        const int gt)
{
    static const unsigned ex[SIZE][N_BASE] =
    {
        { 2, 0, 0, 0},
        { 0, 2, 0, 0},
        { 0, 0, 2, 0},
        { 0, 0, 0, 2},
        { 1, 1, 0, 0},
        { 1, 0, 1, 0},
        { 1, 0, 0, 1},
        { 0, 1, 1, 0},
        { 0, 1, 0, 1},
        { 0, 0, 1, 1}
    };

    return ex[gt][base_id];
}

inline
uint8_t
get_allele(const int gt,
           const unsigned chrom_idx)
{
    static const unsigned ex[SIZE][2] =
    {
        { 0, 0},
        { 1, 1},
        { 2, 2},
        { 3, 3},
        { 0, 1},
        { 0, 2},
        { 0, 3},
        { 1, 2},
        { 1, 3},
        { 2, 3}
    };

    assert(gt<SIZE);
    assert(chrom_idx<2);
    return ex[gt][chrom_idx];
}

/// expect2_bias is a copy of the expect2 function for biased het
/// allele calculations
///
/// 0 is 0.0
/// 1 is het-ratio  (defined by client-code)
/// 2 is (1.-het-ratio)
/// 3 is 1.0
///
/// for consistency, state 1 is always applied to the lower allele
/// value, and state 2 to the higher.
///
inline
unsigned
expect2_bias(const int base_id,
             const int gt)
{
    static const unsigned ex[SIZE][N_BASE] =
    {
        { 3, 0, 0, 0},
        { 0, 3, 0, 0},
        { 0, 0, 3, 0},
        { 0, 0, 0, 3},
        { 1, 2, 0, 0},
        { 1, 0, 2, 0},
        { 1, 0, 0, 2},
        { 0, 1, 2, 0},
        { 0, 1, 0, 2},
        { 0, 0, 1, 2}
    };

    return ex[gt][base_id];
}

inline
const char*
get_vcf_gt(const int gt,
           const int ref_gt)
{
    static const char* gtstr[N_BASE] = { "0/0","0/1","1/1","1/2" };

    static const unsigned ex[SIZE][N_BASE] =
    {
        { 0, 2, 2, 2},
        { 2, 0, 2, 2},
        { 2, 2, 0, 2},
        { 2, 2, 2, 0},
        { 1, 1, 3, 3},
        { 1, 3, 1, 3},
        { 1, 3, 3, 1},
        { 3, 1, 1, 3},
        { 3, 1, 3, 1},
        { 3, 3, 1, 1}
    };

    return gtstr[ex[gt][ref_gt]];
}

/// return diploid gt of two base_ids
// \TODO replace this with a table
inline
unsigned
get_gt_with_alleles(
    const unsigned base1,
    const unsigned base2)
{
    assert(base1 < N_BASE);
    assert(base2 < N_BASE);

    if (base1==base2) return base1;

    for (unsigned gt(N_BASE); gt<SIZE; ++gt)
    {
        if ((DIGT::expect2(base1,gt)==0) ||
            (DIGT::expect2(base2,gt)==0)) continue;
        return gt;
    }

    assert(false && "Can't find diploid genotype");
    return 0;
}

}
