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

/// \file

/// \author Morten Kallberg
///

#pragma once

#include <cmath>
#include <algorithm>
#include <iosfwd>
#include <map>
#include <vector>


struct ranksumObs
{
    ranksumObs() :
        ref(0),
        nonRef(0)
    {}

    unsigned ref;
    unsigned nonRef;
};



// Calculates the Mann-Whitney rank-sum statistic from two populations with
// sparse (and similar) observation spaces
//
class ranksum
{
private:
    std::map<int,ranksumObs> _obsMap;  // observations for ref/alt bases

public:

    // insert an observation for base-call case
    void
    add_observation(
        const bool is_ref,
        const int obs)
    {
        auto& val(_obsMap[obs]);
        if (is_ref)
        {
            val.ref++;
        }
        else
        {
            val.nonRef++;
        }
    }

    //return rank-sum U statistic
    double get_u_stat() const;

    void
    clear()
    {
        _obsMap.clear();
    }

private:
    std::vector<int>
    getSpace() const
    {
        std::vector<int> res;
        for (const auto& val : _obsMap)
        {
            res.push_back(val.first);
        }

        std::sort(res.begin(), res.end());
        return res;
    }

    // returns the count for a given base and observation
    ranksumObs
    get_obs_count(
        const int obs) const
    {
        const auto it = _obsMap.find(obs);
        if (it == _obsMap.end()) return ranksumObs();
        return it->second;
    }

    friend std::ostream& operator<< (std::ostream& out, const ranksum& r);
};

std::ostream& operator<< (std::ostream& out, const ranksum& r);
