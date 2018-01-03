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

#include "sparse_function_1d.hh"

#include <cassert>
#include <cmath>

#ifdef DEBUG_SFUNC
#include <iostream>
#endif



// simple linear function
// construct with 2 points
//
struct linear_function
{

    linear_function(blt_float_t x1,
                    blt_float_t y1,
                    blt_float_t x2,
                    blt_float_t y2)
    {
        assert(std::abs(x2-x1)>0.);
        _s=(y2-y1)/(x2-x1);
        _i=y1-(x1*_s);
    }

    blt_float_t
    get_val(const blt_float_t x) const
    {
        return _i+x*_s;
    }

private:
    blt_float_t _s,_i;
};



inline
bool
in_range(blt_float_t min,
         blt_float_t val,
         blt_float_t max)
{

    return (! ((val<min) || (val>max)));
}



//
blt_float_t
integrate_sparsefunc(const sparse_function& sf,
                     const blt_float_t min_x,
                     const blt_float_t max_x,
                     const blt_float_t min_weight,
                     const blt_float_t max_weight)
{

    assert(min_x < max_x);

    const double range(max_x-min_x);

    assert(range>0);

    assert(min_weight>=0);
    assert(max_weight>=0);
    assert((min_weight+max_weight)>0.);

    // normalize range weights:
    const blt_float_t minw(min_weight/(min_weight+max_weight));
    const blt_float_t maxw(max_weight/(min_weight+max_weight));

    const linear_function wf(min_x,minw,max_x,maxw);

    blt_float_t sum(0.);

    typedef sparse_function::const_iterator sit;

    sit i(sf.begin());
    const sit i_end(sf.end());

    bool is_last_set(false);
    double last_x(0);
    double last_y(0);

    // use values outside of range for one step only -- don't use an
    // interval that's fully off range.
    //
    // the reasoning behind this is that we don't want to lose 0 and 1
    // values due to any type of float rounding issue:
    //

    // write first version without handling boundaries correctly

    // handled covered boundaries:

    // handle uncovered boundaries:

    //
    for (; i!=i_end; ++i)
    {
        const blt_float_t& x(i->first);
        const blt_float_t& y(i->second);

        const bool is_in_range(in_range(min_x,x,max_x));
        bool is_last_in_range(false);

        if (! is_last_set)
        {
            if (is_in_range)
            {
                // implicit start segment:
                const blt_float_t segrange(x-min_x);
                sum += (y*(wf.get_val(min_x)+wf.get_val(x))) * (segrange/range);
            }

        }
        else
        {
            is_last_in_range=(in_range(min_x,last_x,max_x));

            if (is_in_range &&
                is_last_in_range)
            {
                // entire segment is enclosed in target range:

                const blt_float_t segrange(x-last_x);
                sum += ((last_y*wf.get_val(last_x))+(y*wf.get_val(x))) * (segrange/range);

            }
            else if (is_in_range)
            {
                // segment enters target range:

                const blt_float_t segrange(x-min_x);
                const blt_float_t func_segrange(x-last_x);

                if (func_segrange>0.)
                {
                    const blt_float_t frac(segrange/func_segrange);
                    assert(frac<=1.);
                    // simple linear interpolate:
                    const blt_float_t min_y(y*(1.-frac)+last_y*(frac));
                    sum += ((min_y*wf.get_val(min_x))+(y*wf.get_val(x))) * (segrange/range);
                }
            }
            else if (is_last_in_range)
            {
                // segment exits target range:

                const blt_float_t segrange(max_x-last_x);
                const blt_float_t func_segrange(x-last_x);

                if (func_segrange>0.)
                {
                    const blt_float_t frac(segrange/func_segrange);
                    assert(frac<=1.);
                    // simple linear interpolate:
                    const blt_float_t max_y(last_y*(1.-frac)+y*(frac));
                    sum += ((max_y*wf.get_val(max_x))+(last_y*wf.get_val(last_x))) * (segrange/range);
                }
            }
        }

        last_x=x;
        last_y=y;
        is_last_set=true;

        if (is_last_in_range && (! is_in_range)) break;
    }

    assert(is_last_set);

    if (in_range(min_x,last_x,max_x))
    {
        // implicit exit segment:
        const blt_float_t segrange(max_x-last_x);
        sum += (last_y*(wf.get_val(max_x)+wf.get_val(last_x))) * (segrange/range);
    }

    return sum;
}



// sparse function is scaled out of log representation, then integrated.
//
blt_float_t
integrate_ln_sparsefunc(const sparse_function& sf,
                        const blt_float_t min_x,
                        const blt_float_t max_x,
                        const blt_float_t min_weight,
                        const blt_float_t max_weight)
{

    typedef sparse_function::const_iterator sit;

    const sit i_start(sf.begin()),i_end(sf.end());

    bool is_max(false);
    blt_float_t max(0);
    for (sit i(i_start); i!=i_end; ++i)
    {
        if ((! is_max) || (i->second > max))
        {
            max=i->second;
            is_max=true;
        }
    }

#ifdef DEBUG_SFUNC
    std::cerr << "SFUNC: max: " << max << "\n";
#endif

    sparse_function scaled_sf;
    for (sit i(i_start); i!=i_end; ++i)
    {
        scaled_sf.insert(i->first,std::exp(i->second-max));
#ifdef DEBUG_SFUNC
        std::cerr << "SFUNC: sf: x,y: " << i->first << " " << i->second << "\n";
        std::cerr << "SFUNC: ssf: x,y: " << i->first << " " << std::exp(i->second-max)  << "\n";
#endif
    }

    const blt_float_t val(integrate_sparsefunc(scaled_sf,min_x,max_x,min_weight,max_weight));

#ifdef DEBUG_SFUNC
    std::cerr << "WAG: val " << val << "\n";
#endif

    return std::log(val)+max;
}
