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
/// \author Morten Kallberg
///

#include "blt_util/fastRanksum.hh"

#include <cmath>
#include <iostream>



static
double
getMannWhitneyZScore(
    const int n1,
    const int n2,
    const double w1)
{
    const double mean = n1*(n1+n2+1)/2.0;
    const double var  = std::sqrt(n2*mean/6.0);
    static const double epsilon(0.0001);
    if (std::abs(var)<epsilon)
    {
        return 0.0;
    }
    return (w1-mean)/var;
}



double
fastRanksum::
getExpectedCategory2Value() const
{
    double R2 = 0;
    unsigned N2 = 0;

    for (unsigned i=0; i<_obs.size(); i++)
    {
        if (_obs[i].empty()) continue;
        const unsigned c2 = _obs[i].c2;
        R2 += i*c2;
        N2 += c2;
    }

    if (N2==0) return 0;
    return R2/N2;
}



double
fastRanksum::
get_z_stat() const
{
    double R1 = 0;
    double R2 = 0;
    int N1 = 0;
    int N2 = 0;

    for (const auto& robs : _obs)
    {
        if (robs.empty()) continue;
        const int obs1 = robs.c1;
        const int obs2 = robs.c2;
        double rank_weight = (2*(N1+N2+1) + (obs1+obs2)-1)/2.0;
        R1 += rank_weight*obs1;
        R2 += rank_weight*obs2;
        N1 += obs1;
        N2 += obs2;
    }

#if 1
    // this performs the equivalent of the z-score computation on U, even though U, mean(U), var(U)
    // are not directly enumerated
    // TODO: is there a reference for this transformation?
    if (R1>R2)
    {
        return getMannWhitneyZScore(N2, N1, R2);
    }
    else
    {
        return getMannWhitneyZScore(N1, N2, R1);
    }
#else
    // for verification purposes, you can actually enumerate U and U stats to get the ?always same? answer here:
    const double U1 = R1 - (N1*(N1+1))/2.;
    const double U2 = R2 - (N2*(N2+1))/2.;
    const double U = std::min(U1,U2);

    const double mean = N1*N2/2.0;
    const double var  = std::sqrt(mean*(N1+N2+1)/6.0);

    static const double epsilon(0.0001);
    if (std::abs(var)<epsilon)
    {
        return 0.0;
    }
    return (U-mean)/var;
#endif
}
