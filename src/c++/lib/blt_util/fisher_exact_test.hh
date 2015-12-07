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
/// \author Mitch Bekritsky
/// \author Peter Krusche
///

#pragma once

enum class FISHER_EXACT : unsigned int {
    LESS      = 1,
    GREATER   = 2,
    TWOTAILED = 3,
};

/**
 * Compute the fisher exact test p-value for a 2x2 contingency matrix:
 *            group     P      Q
 * condition
 *        C             a      b
 *        D             c      d
 *
 * We do a two-tailed test.
 */
extern
double
fisher_exact_test_pval_2x2(
    const unsigned a,
    const unsigned b,
    const unsigned c,
    const unsigned d,
    const FISHER_EXACT type=FISHER_EXACT::TWOTAILED);
