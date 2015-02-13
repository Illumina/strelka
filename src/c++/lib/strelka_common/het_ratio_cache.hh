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

/// \author Chris Saunders
///

#pragma once

#include "blt_util/blt_types.hh"

#include <array>
#include <vector>


/// A simple static sized array with deep copy semantics:
///
template <unsigned NVAL>
struct cache_val
{
    std::array<blt_float_t,NVAL> val;
};


/// a special-case memoizer for lhood computation intermediates
///
/// Effectively func(qscore,hetratio) -> array[NVAL] is being memoized
///
/// the special case logic here is that we know both qscore and het ratio are
/// bound to ranges on [0,small int] so we can lookup in a vector instead
/// of some type of tree/hash
///
/// The interface is totally non-general and pushes all of the work onto the client (sorry):
///
/// get_val() returns a tuple, the first value of which is
/// a bool indicating whether the returned data structure has already
/// been called, and thus is (presumably) cached. Note that the client is
/// responsible for setting any values into the returned array for
/// caching.
///
///
/// this value caching didn't do much for the grid model -- better to
/// leave it out for now... (how up to date is this comment?)
//
// TODOs:
//  - how can we take advantage of this with a more maintainable interface:
//      - func(x), memo(func(x)...
//  - another bounded memoizer in the code is used for qscores -- is there a generalized interface:
//      - small_index_memoizer?
//  - starling could use this if we get rid of the dependent eprob business
//
//
template <unsigned NVAL>
struct het_ratio_cache
{
    het_ratio_cache()
        : _is_cached(MAX_QSCORE* MAX_INDEX,false)
        , _cache(MAX_QSCORE* MAX_INDEX)
    {}

    std::pair<bool,cache_val<NVAL>*>
    get_val(const unsigned qscore,
            const unsigned ratio_index)
    {
        if (qscore>=MAX_QSCORE ||
            ratio_index>=MAX_INDEX)
        {
            return std::make_pair(false,&_any_val);
        }

        const unsigned index(ratio_index + qscore*MAX_INDEX);
        if (_is_cached[index])
        {
            return std::make_pair(true,&(_cache[index]));
        }
        else
        {
            _is_cached[index] = true;
            return std::make_pair(false,&(_cache[index]));
        }
    }

private:
    enum contanst { MAX_QSCORE = 64, MAX_INDEX = 12 };

    typedef cache_val<NVAL> cache_val_n;

    cache_val_n _any_val; // return this if a request is outside of cached range
    std::vector<bool> _is_cached;
    std::vector<cache_val_n> _cache;
};
