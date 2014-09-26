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
#include <iosfwd>
#include <iostream>

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
//        if (!isA){
//            this->total_alt+=obs;
//            this->total_alt_count++;
//            std::cerr << "obs: " << obs << std::endl;
//            std::cerr << "sum: " <<this->total_alt << std::endl;
//            std::cerr << "count: " << this->total_alt_count << std::endl;
//        }
    }

    //return rank-sum U statistic
    double get_u_stat() const;

    //
    double get_u_stat_uniform() const;

    // return average cycle of the alternate allele
    double get_raw_score() const;

    double get_avg_alt() const;

    void
    clear()
    {
        _obs.resize(0);
//        this->total_alt = 0;
//        this->total_alt_count = 0;
    }


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
            if (isA) A++; // case
            else
            {
                B++;
            }
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
//    unsigned total_alt=0;
//    unsigned total_alt_count=0;
};
