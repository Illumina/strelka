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

/// \author Chris Saunders
/// \author Mitch Bekritsky
/// \author Peter Krusche
///

#pragma once

enum class FISHER_EXACT : unsigned int
{
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
double
fisher_exact_test_pval_2x2(
    const unsigned a,
    const unsigned b,
    const unsigned c,
    const unsigned d,
    const FISHER_EXACT type=FISHER_EXACT::TWOTAILED);
