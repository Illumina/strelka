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
        return std::move(res);
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
