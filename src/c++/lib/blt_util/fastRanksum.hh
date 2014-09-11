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

///
/// \author Morten Kallberg
///

#pragma once

#include <vector>


// Calculates the Mann-Whitney rank-sum statistic from two populations with
// sparse (and similar) observation spaces
//
// this is a variation on the original ranksum class which only accepts small unsigned
// observations for better performance.
//
struct fastRanksum
{
    // insert an (A/B) observation for base-call case
    void
    add_observation(
        const bool isA,
        const unsigned obs)
    {
        if (obs >= _obs.size()) _obs.resize(obs+16);
        _obs[obs].inc(isA);
    }

    //return rank-sum U statistic
    double get_u_stat() const;

    void
    clear()
    {
        _obs.resize(0);
    }

    double get_raw_score();

private:

    struct ranksumObs
    {
        ranksumObs() :
            A(0),
            B(0)
        {}

        void
        inc(const bool isA)
        {
            if (isA) A++;
            else     B++;
        }

        bool
        empty() const
        {
            return ((A+B)==0);
        }

        unsigned A;
        unsigned B;
    };

    std::vector<ranksumObs> _obs; // observations for ref/alt bases
};
