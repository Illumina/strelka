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
/// \author Mitch Bekritsky
/// \author Chris Saunders
///

#include "min_count_binom_gte_cache.hh"

#include "blt_util/binomial_test.hh"

#include "boost/math/special_functions/gamma.hpp"

#include <cassert>



/// TODO this is bulky -- isn't there an easier way to have an array size constant in C++11?
constexpr unsigned min_count_binom_gte_cache::_maxCacheK;



min_count_binom_gte_cache::
min_count_binom_gte_cache(
    const double alpha)
    : _alpha(alpha)
{
    assert((alpha >= 0.) && (alpha <= 1.));

    // caching scheme relies on conditions where binom is approximated by Poisson, so
    // X ~ B(n,p) is approximated by X ~ Poisson(\lambda), where \lambda is n*p. This
    // approximation allows us (with fixed \alpha), to cache in a compact 1D space.
    //
    // For a given \alpha we want to know for what value of X = k we exceed \alpha on the
    // Poisson CDF: F(k,\lambda). Because we are caching responses here for unknown
    // values of \lambda, we want the inverse of the CDF. By first expressing
    // the CDF, F(k, \lambda), as the equivalent Q(k+1,\lambda), where Q is the (upper)
    // incomplete gamma function, we can then use the inverse gamma function to
    // provide \lambda = Q^{-1}(k+1,\alpha)
    //
    // In practice boost supplies Q^{-1}() via boost::math::gamma_q_inv
    //
    // A helpful reference to work this out comes from stackexchange:
    // http://math.stackexchange.com/questions/68258/poisson-cdf-and-solving-for-lambda
    //
    // Note that in practice the inverse upper incomplete gamma gives 1-x, where x is the
    // expected answer based on the writeup and external links above, so the implementation
    // below uses the inverse lower incomplete gamma, provided by boost as boost::math::gamma_p_inv
    //
    /// \TODO explain the unexpected complement of the expected value when using
    ///    gamma_q_inv and fix the documentation above to be consistent
    //
    for (unsigned kminus1(0); kminus1<_maxCacheK; ++kminus1)
    {
        _papprox[kminus1] = boost::math::gamma_p_inv(kminus1+1,alpha);
    }
}


static
bool
isPoissonApprox(
    const double p)
{
    // Test whether Poisson approx to binomial is reasonable for this (n,p) by checking
    // that the Poisson/bionomial variance ratio is close to 1
    //
    // ratio of Pois to Binom variance is var_ratio = n*p / (n*p*(1-p)) -> 1/(1-p)
    // comp_var_ratio = 1 - var_ratio -> 1 - 1/(1-p) -> -p/(1-p)
    const double abs_comp_var_ratio(p/(1-p));
    return (abs_comp_var_ratio <= 0.01);
}



unsigned
min_count_binom_gte_cache::
getCount(
    const unsigned n,
    const double p) const
{
    assert(n>0);
    assert(p>=0. && p<=1.);

    if (p<=0.) return 1.;

    if (isPoissonApprox(p))
    {
        const auto iter = std::lower_bound(_papprox.begin(),_papprox.end(),n*p);
        if (iter != _papprox.end())
        {
            return (1+(iter-_papprox.begin()));
        }
    }
    return min_count_binomial_gte_exact(_alpha,p,n);
}



// Worksheet to capture the transformation here:
//
// n_success >= getCount(n,p)
// n_success >= index+1, where pa[index] > np
// n_success > index, where pa[index] > np
// pa[n_success-1] > np
//
// usually faster than calling getCount b/c we don't need
// to search on the papprox array
//
bool
min_count_binom_gte_cache::
isRejectNull(
    const unsigned n,
    const double p,
    const unsigned n_success) const
{
    assert(n>0);
    assert(p>=0. && p<=1.);
    assert(n_success<=n);

    if (p<=0.) return (n_success>0);
    if (n_success==0) return false;

    if (isPoissonApprox(p))
    {
        const double np(n*p);
        if (np < _papprox.back())
        {
            const unsigned tn(std::min(n_success,_maxCacheK)-1);
            return (np < _papprox[tn]);
        }
    }

    /// TODO I'm sure there's a faster test encapsulating the line below into a single isRejectNull predicate
    return (n_success>=min_count_binomial_gte_exact(_alpha,p,n));
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
