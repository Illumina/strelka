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
///

#pragma once

#include <array>


/// \brief special caching variant of the min_count_binomial_gte_exact function
///
/// The min_count_binomial_gte_exact function returns the minimum number of
/// successes to reject the null hypothesis with a p-value of at most alpha
/// for a given error rate (p) and number of trials (n).
///
/// This object is a special caching variant of the min_count_binomial_gte_exact
/// function. For caching to be effective you must (1) supply a single value of
/// alpha beforehand and (2) The (n,p) range of your binomial calls should meet
/// Poisson approximation conditions. The optimized implementation takes advantage
/// of the Poisson approximation to the binomial to efficiently cache the 2D (n,p)
/// function space to a more compact 1D (n*p) space. Calls to (n,p) which fall outside
/// the range of Poisson approximation will be routed to the standard
/// min_count_binomial_gte_exact function, so results for this case will be accurate
/// but have unoptimized compute cost.
///
struct min_count_binom_gte_cache
{
    explicit
    min_count_binom_gte_cache(const double alpha);

    /// matches min_count_binomial_gte_exact function
    unsigned
    getCount(
        const unsigned n,
        const double p) const;

    /// \returns equiv of (n_success >= getCount(n,p))
    ///
    /// this version is faster than calling getCount() and doing
    /// your own comparison when approx coniditions hold
    bool
    isRejectNull(
        const unsigned n,
        const double p,
        const unsigned n_success) const;

private:
    unsigned
    get_min_count_approx(const double np) const;

    const double _alpha;

    static constexpr unsigned _maxCacheK = 18;
    std::array<double,_maxCacheK> _papprox;
};
