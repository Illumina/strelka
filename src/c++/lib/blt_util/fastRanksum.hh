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

#pragma once

#include <vector>
#include <iosfwd>


/// Calculates the Mann-Whitney rank-sum statistic from two populations with
/// sparse (and similar) observation spaces
///
/// this is a variation on the original ranksum class which only accepts small unsigned
/// observations for better performance.
///
struct fastRanksum
{
    /// insert an observation, indicating membership in category 1 or 2
    void
    add_observation(
        const bool isCategory1,
        const unsigned obs)
    {
        if (obs >= _obs.size()) _obs.resize(obs+16);
        _obs[obs].inc(isCategory1);
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

    /// \return z-score of the Mann-Whitney U statistic
    double get_z_stat() const;

    /// \return average value of category 2
    double getExpectedCategory2Value() const;

    void
    clear()
    {
        _obs.resize(0);
    }


private:

    struct ranksumObs
    {
        ranksumObs() :
            c1(0),
            c2(0)
        {}

        void
        inc(const bool isCategory1)
        {
            if (isCategory1) c1++;
            else
            {
                c2++;
            }
        }

        void
        merge(
            const ranksumObs& rhs,
            const bool isFlipObservationCategories)
        {
            if (isFlipObservationCategories)
            {
                c1 += rhs.c1;
                c2 += rhs.c2;
            }
            else
            {
                c1 += rhs.c2;
                c2 += rhs.c1;
            }
        }

        bool
        empty() const
        {
            return ((c1+c2)==0);
        }

        unsigned c1;
        unsigned c2;
    };

    std::vector<ranksumObs> _obs;
};
