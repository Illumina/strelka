// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
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
#include <cmath>
#include <algorithm>
#include <iosfwd>
#include <map>
#include <vector>


#pragma once


// Calculates the Mann-Whitney rank-sum statistic from two populations with
// sparse (and similar) observation spaces
//
class ranksum
{
private:
    char ref_base;

public:
    int N1;
    int N2;
    double R1;
    double R2;
    std::map<int,int> l1;  // observations for reference base
    std::map<int,int> l2;  // observations for alternate base
    std::map<int,int> space;// full key space
    //constructor

    ranksum() {
        set_ref_base('N');
    }

    ranksum(char base) {
        set_ref_base(base);
//    	N1(0);
//		N2(0);
//		U1(0.0);
//		U2(0.0);
//		R1(0.0);
//		R2(0.0);
    }

    // insert an observation for base-call case
    void add_observation(bool is_ref, int obs);

    //return rank-sum U statistic
    double get_u_stat();

    void set_ref_base(char base) {
        ref_base = base;
    }

    std::vector<int> getSpace() {
        std::vector<int> res;
        for (std::map<int, int>::iterator it = space.begin(); it != space.end(); ++it) {
            res.push_back((*it).first);
        }

        std::sort(res.begin(), res.end());
        return res;
    }

    // returns the count for a given base and observation
    int get_obs_count(bool is_ref, int obs) {

        int res = 0;
        std::map<int,int>::iterator it;

        //check for observation in reference map
        if (is_ref) {
            it = l1.find(obs);
            if (it != l1.end())
                res = it->second;
        }
        else {
            it = l2.find(obs);
            if (it != l2.end())
                res = it->second;
        }
        return res;
    }

    void clear() {
        l1.clear();
        l2.clear();
        space.clear();

    }

    char get_refbase() {
        return ref_base;
    }
};

std::ostream& operator<< (std::ostream& out, ranksum& r);
