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

/// \author Chris Saunders
///
#ifndef __SPARSE_FUNCTION_1D_HH
#define __SPARSE_FUNCTION_1D_HH

#include "blt_util/blt_types.hh"

#include <map>



// this object holds values of a function sampled at non-regular
// intervals, it allows the following operations:
//
// (1) insert new value: f(x) = y
//
// (2) iterate through all existing values sorted on x
//
struct sparse_function
{

    typedef std::map<blt_float_t,blt_float_t> fmap_t;
    typedef fmap_t::const_iterator const_iterator;

    const_iterator
    begin() const
    {
        return _fmap.begin();
    }

    const_iterator
    end() const
    {
        return _fmap.end();
    }

    // we don't check for repeated x:
    //
    void
    insert(const blt_float_t x,
           const blt_float_t y)
    {
        _fmap.insert(std::make_pair(x,y));
    }

private:
    fmap_t _fmap;
};



// Approximate integration of f(x)P(x), f(x) is sparse, over the
// specified x-range assuming P(x) is as simple linear slope over
// the range defined by min and max weight. If min and max weight are
// set to 1 a uniform distribution is used.
//
// method assumes a linear connection between sampled points in the
// function, including points outside of the min,max range. when
// sampling does not extend to the end of either range the last
// sampled value is treated as constant for the remainder of the range
//
blt_float_t
integrate_sparsefunc(const sparse_function& sf,
                     const blt_float_t min_x,
                     const blt_float_t max_x,
                     const blt_float_t min_weight,
                     const blt_float_t max_weight);


// Approximate integration of ln(f(x))P(x), ln(f(x)) is sparse, over
// the specified x-range assuming P(x) is as simple linear slope over
// the range defined by min and max weight. If min and max weight are
// set to 1 a uniform distribution is used.
//
// sparse function is scaled out of log representation, then integrated.
//
blt_float_t
integrate_ln_sparsefunc(const sparse_function& sf,
                        const blt_float_t min_x,
                        const blt_float_t max_x,
                        const blt_float_t min_weight,
                        const blt_float_t max_weight);

#endif
