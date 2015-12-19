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
///

#pragma once

#include <vector>


/// special caching variant of min_count_binomial_gte_exact
///
/// min_count_binomial_gte_exact returns the minimum number of successes to
/// reject the null hypothesis with a p-value of at most alpha for a given
/// error rate and number of trials.
///
/// this is a special caching variant of that function which, for a given
/// value of alpha takes advantage of the Poisson approximation to the binomial
/// to efficiently cache the 2D (n,p) function space to the 1D (np)
/// space. For this reason the caching used by this method is only effective
/// under Poisson approximation conditions (n 'large'/ p 'small').
///
struct min_count_binom_gte_cache
{
    explicit
    min_count_binom_gte_cache(const double alpha);

    unsigned
    getCount(
        const unsigned n,
        const double p) const;

private:
    unsigned
    get_min_count_approx(const double np) const;

    const double _alpha;
    std::vector<double> _papprox;
};
