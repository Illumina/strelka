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
/// \author Chris Saunders
///

#include <cassert>


struct fval
{
    fval(const blt_float_t x_init = 0.,
         const blt_float_t val_init = 0.)
        : x(x_init)
        , val(val_init) {}

    blt_float_t x;
    blt_float_t val;
};



// Helper function to create new fval object from x value, while
// recording the sampled function value. Side-effects are too heavy to
// make a ctor out of this:
template <typename Func>
fval
make_fval(const blt_float_t x,
          Func& f)
{
    const blt_float_t val(f.calc_and_store_val(x));
    return fval(x,val);
}



template <typename Func>
void
sample_uniform_range(const blt_float_t min_x,
                     const blt_float_t max_x,
                     Func& f)
{

    static const blt_float_t min_range(0.0001);
    static const unsigned max_iter(8);

    assert(min_x < max_x);

    fval low(make_fval(min_x,f));
    fval high(make_fval(max_x,f));

    for (unsigned i(0); i<max_iter; ++i)
    {

        const blt_float_t range(high.x-low.x);
        const blt_float_t next_range(range*0.5);

        assert(next_range>=0);
        if (next_range<min_range) return;

        const fval mid(make_fval(low.x+next_range,f));

        if (high.val>low.val)
        {
            if ((high.x+next_range)<max_x)
            {
                const fval mid2(make_fval(high.x+next_range,f));
                if (mid2.val>mid.val)
                {
                    low=high;
                    high=mid2;
                }
                else
                {
                    low=mid;
                }
            }
            else
            {
                low=mid;
            }
        }
        else
        {
            if ((low.x-next_range)>min_x)
            {
                const fval mid2(make_fval(low.x-next_range,f));
                if (mid2.val>=mid.val)
                {
                    high=low;
                    low=mid2;
                }
                else
                {
                    high=mid;
                }
            }
            else
            {
                high=mid;
            }
        }
    }
}
