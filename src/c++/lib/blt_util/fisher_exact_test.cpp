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

#include "fisher_exact_test.hh"

#include <boost/math/distributions/hypergeometric.hpp>
#include <algorithm>

using namespace boost::math;
using namespace std;

double
fisher_exact_test_pval_2x2(
    const unsigned a,
    const unsigned b,
    const unsigned c,
    const unsigned d,
    const FISHER_EXACT type)
{
    /**
     *             defective  good
     *  sampled    a          b
     *  remaining  c          d
     *
     *  We want to know the probability of having drawn a defective and b good
     *  items without replacement.
     */

    // http://mathworld.wolfram.com/FishersExactTest.html
    const unsigned N = a + b + c + d;  // total number of elements in urn
    const unsigned r = a + c;          // total number of defective elements
    const unsigned n = a + b;          // total number of elements we sampled
    const unsigned k = a;              // number of defective elements we got in our sample

    unsigned max_for_k = min(r, n);   // maximum number of defective elements we can draw in a sample of size n
    unsigned min_for_k = (unsigned )max(0, int(r + n - N)); // once we're out of good elements, the remaining ones
    // are defective.

    if (type == FISHER_EXACT::LESS)
    {
        max_for_k = a;
    }

    if (type == FISHER_EXACT::GREATER)
    {
        min_for_k = a;
    }

    hypergeometric_distribution<> hgd(r, n, N);

    // for the two-tailed test, sum over both tails where p <= p_cutoff
    // for the one-tailed tests, we have excluded tails of the distribution
    // by changing the limits of summation above
    //

    // Per the implementation of 'fisher.test' in R, note that the cutoff needs an epsilon term to account
    // for minor floating point discrepancies across platforms. In particular, this has been found to be
    // important for long double math differences between linux and windows/WSL. WSL details of this
    // issue are described here:
    //
    // https://github.com/Microsoft/WSL/issues/830
    //
    // An example input triggering the issue on window is an (a,b,c,d,type) input of (54,107,0,2,twotailed),
    // where the expected result is 1
    //
    static const double relativeError(1. + 1.e-7);
    const double cutoff(relativeError * ((type == FISHER_EXACT::TWOTAILED) ? pdf(hgd, k) : 1.0) );
    double p = 0;
    for (unsigned q = min_for_k; q <= max_for_k; ++q)
    {
        const double _p = pdf(hgd, q);
        if (_p <= cutoff) p += _p;
    }
    return p;
}
