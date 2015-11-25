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

/// \file

/// \author Morten Kallberg
///

#include "blt_util/fastRanksum.hh"

#include <cmath>
#include <iostream>



static
double
get_z_score(const int n1, const int n2, const double w1)
{
    const double mean = n1*(n1+n2+1)/2.0;
    const double var  = std::sqrt(n2*mean/6.0);
    if (static_cast<int>(var)==0)
    {
        return 0.0;
    }
    return (w1-mean)/var;
}

double
fastRanksum::
get_raw_score() const
{
    double R2 = 0;
    int N2 = 0;

    for (unsigned i=0; i<_obs.size(); i++)
    {
        if (_obs[i].empty()) continue;
        const int obs2 = _obs[i].B;
        R2 += i*obs2;
        N2 += obs2;
    }

    //return the z-score for the smaller of U1 and U2 assuming a gaussian approx.
    if (N2>0)
        return 1.0*R2/N2;
    return 0;
    //    return get_avg_alt();
}

// return the U statistic
double
fastRanksum::
get_u_stat() const
{
    double R1 = 0;
    double R2 = 0;
    int N1 = 0;
    int N2 = 0;

    for (const auto& robs : _obs)
    {
        if (robs.empty()) continue;
        const int obs1 = robs.A;
        const int obs2 = robs.B;
        double rank_weight = (2*(N1+N2+1) + (obs1+obs2)-1)/2.0;
        R1 += rank_weight*obs1;
        R2 += rank_weight*obs2;
        N1 += obs1;
        N2 += obs2;
    }
    //return the z-score for the smaller of U1 and U2 assuming a gaussian approx.
    if (R1>R2)
    {
        return get_z_score(N2,N1,R2);
    }
    else
    {
        return get_z_score(N1,N2,R1);
    }
}

double
fastRanksum::
get_u_stat_uniform() const
{
    double R1 = 0;
    double R2 = 0;
    int N1 = 0;
    int N2 = 0;
    //loop over all observations
    for (unsigned i=0; i<_obs.size(); i++)
    {
        if (_obs[i].empty()) continue;
        const int obs1 = _obs[i].A;
        const int obs2 = _obs[i].B;
//        double rank_weight = (2*(N1+N2+1) + (obs1+obs2)-1)/2.0;
        R1 += i*obs1;
        R2 += i*obs2;
        N1 += obs1;
        N2 += obs2;
    }
    //return the average rank weight of case with most observations
    if (N1>N2)
        return R1/N1;
    return R2/N2;
}
