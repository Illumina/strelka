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

/// \file

/// \author Morten Kallberg
///

#include "blt_util/fastRanksum.hh"

#include <cmath>

#include "blt_util/log.hh"


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
get_raw_score() const{
    double R1 = 0;
    double R2 = 0;
    int N1 = 0;
    int N2 = 0;

    //loop over all observations
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

    //return the average rank weight of case with most observations
    if (N1>N2)
        return R1/N1;
    return R2/N2;
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

    //loop over all observations
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
