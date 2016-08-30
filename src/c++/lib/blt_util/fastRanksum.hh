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
/// \author Morten Kallberg
///

#pragma once

#include <vector>
#include <iosfwd>


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

    void
    merge(
        const fastRanksum& rhs,
        const bool isFlipObservationCategories = false)
    {
        if (rhs._obs.size() > _obs.size()) _obs.resize(rhs._obs.size());
        const unsigned obsCount(_obs.size());
        for (unsigned obsIndex(0); obsIndex < obsCount; ++obsIndex)
        {
            _obs[obsIndex].merge(rhs._obs[obsIndex], isFlipObservationCategories);
        }
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

        void
        merge(
            const ranksumObs& rhs,
            const bool isFlipObservationCategories)
        {
            if (isFlipObservationCategories)
            {
                A += rhs.A;
                B += rhs.B;
            }
            else
            {
                A += rhs.B;
                B += rhs.A;
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
