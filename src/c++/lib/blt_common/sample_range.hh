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

#pragma once


// Interface for sample function object, sampled function should be
// able to inherit from this interface.
//
#if 0
struct sample_func_iface
{
    virtual
    blt_float_t
    calc_and_store_val(const blt_float_t x) = 0;
};
#endif


// A cheap deterministic pseudo-importance sampler
//
// Object func must meet requirements described above for
// sample_func_iface
//
template <typename Func>
void
sample_uniform_range(const blt_float_t min_x,
                     const blt_float_t max_x,
                     Func& f);


#include "sample_range_impl.hh"
