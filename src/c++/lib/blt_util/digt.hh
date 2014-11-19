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

///
/// \author Chris Saunders
///

#pragma once

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
    switch (idx)
    {
    case AA:
    case CC:
    case GG:
    case TT:
        return false;
    default:
        return true;
    }
}

inline
char
hom_base(const unsigned idx)
{
    switch (idx)
    {
    case AA:
        return 'A';
    case CC:
        return 'C';
    case GG:
        return 'G';
    case TT:
        return 'T';
    default:
        return 'X';
    }
}

/// provide the alelle frequency for base_id given a diploid genotype gt
inline
double
expect(const int base_id,
       const int gt)
{
    static const unsigned N_BASE(4);

    static const double ex[SIZE][N_BASE] = {{ 1.0, 0.0, 0.0, 0.0},
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
    static const unsigned N_BASE(4);

    static const unsigned ex[SIZE][N_BASE] = {{ 2, 0, 0, 0},
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
    static const unsigned N_BASE(4);

    static const unsigned ex[SIZE][N_BASE] = {{ 3, 0, 0, 0},
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
    static const unsigned N_BASE(4);

    static const char* gtstr[N_BASE] = { "0/0","0/1","1/1","1/2" };

    static const unsigned ex[SIZE][N_BASE] = {{ 0, 2, 2, 2},
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
    static const unsigned N_BASE(4);
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
