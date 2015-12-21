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
/// \author Mitch Bekritsky
/// \author Chris Saunders
///

#include "min_count_binom_gte_cache.hh"

#include "blt_util/binomial_test.hh"

#include <boost/math/special_functions/gamma.hpp>

#include <cassert>



min_count_binom_gte_cache::
min_count_binom_gte_cache(
    const double alpha)
  : _alpha(alpha)
{
    // caching scheme relies on conditions where binom is approximated by Poisson, so
    // X ~ B(n,p) -> X ~ Poisson(\lambda), where \lambda is n*p
    //
    // For a given alpha we want to know for what value of k we exceed alpha on the
    // Poisson cdf: F(k,\lambda). Because we are caching responses here in response to
    // unknown values of \lambda, we want the inverse of F(k, \lambda). By first tranforming to
    //
    // F(k, \lambda) -> P(k+1,\lambda), where P is the (lower) incomplete gamma function, we then
    // rely on inversions of the the gamma function to provide \lambda = P(k+1,alpha)^{-1}
    //
    // In practice boost supplies P()^{-1} via boost::math::gamma_p_inv
    //
    // Besides wikipedia, a critical reference to work this out comes from stackexchange here:
    //
    // http://math.stackexchange.com/questions/68258/poisson-cdf-and-solving-for-lambda
    //
    static const unsigned maxCacheK(18);
    for (unsigned kminus1(0);kminus1<maxCacheK;++kminus1)
    {
        _papprox.push_back(boost::math::gamma_p_inv(kminus1+2,alpha));
    }
}


unsigned
min_count_binom_gte_cache::
get_min_count_approx(const double np) const
{
    // this should be faster than a binary interval search since
    // we almost always expect n * p to be on the low side
    // (i.e. np < 0.01)
    const unsigned pSize(_papprox.size());
    for(unsigned kminus1 = 0; kminus1 < pSize; ++kminus1)
    {
        if(np > _papprox[kminus1]) continue;
        return (kminus1+1);
    }
    return 0;
}



unsigned
min_count_binom_gte_cache::
getCount(
    const unsigned n,
    const double p) const
{
    assert(n>0);
    assert(p>=0.);

    if (p<=0.) return 1.;

    /// TODO: this is minimal disaster-proofing, insert a real test for validity of poisson approx to binomial.
    ///    also low n is common for the application we have in mind, how to improve this case?
    double mv_ratio = (n * p) / (n * (1 - p));
    if (abs(1 - mv_ratio) <= 0.01)
    {
        unsigned val = get_min_count_approx(n*p);
        if (val != 0) return val;
    }
    return min_count_binomial_gte_exact(_alpha,p,n);
}


// original docs form Mitch:

// we can leverage the fact the n * p ~ n * p * (1 - p)
// for our current indel error rates and use a Poisson approximation
// to the binomial distribution.   This will not work when p is large
// relative to n (i.e. high error rates--in a typical 30x dataset, error
// rates would have to be around 10%, we're currently observing error rates
// around 5e-6 - 2e-5.  This lookup table provides the minimum coverage
// given an error rate p and coverage n when for a p-value of 1e-9.
// We could calculate this dynamically in a small piece of code while we
// do refCoverage, etc. for an arbitrary p-value, but this should do for now.

/*
* Generated from the following code in R
*   require(plyr)

*   min_success_enum <- data.frame(min_success = qpois(1e-9, lambda = seq(1e-6, 3, 1e-6),
*                                                      lower.tail = FALSE),
*                                  np = seq(1e-6, 3, 1e-6))
*   min_success_enum$min_success[min_success_enum$min_success < 2] <- 2

*   min_success_ranges <- ddply(min_success_enum, .(min_success), summarize,
*                               min_success = unique(min_success),
*                               min_np = min(np),
*                               max_np = max(np))
*   min_success_ranges$min_np <- min_success_ranges$min_np - 1e-6

*   cat(apply(min_success_ranges, 1, function(x) sprintf("{%6d, %.6f, %.6f},",
*                                             x["min_success"],
*                                             x["min_np"], x["max_np"])),
*       sep = "\n")

*/
