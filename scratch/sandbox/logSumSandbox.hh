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

/// \file
/// \brief Alternate log sum implementations
///
/// This file contains various working log sum prototypes
///

#pragma once

#include "blt_util/math_util.hh"


/// Variadic version of getLogSum... this is slightly less fast (~20%) than the version in logSumUtil coded for a
/// specific number of arguments
//
template<typename T0, typename... Ts>
typename std::common_type<T0, Ts...>::type
getLogSum(
    T0 t,
    Ts... vs)
{
    typedef typename std::common_type<T0, Ts...>::type value_type;
    static_assert(std::is_floating_point<value_type>::value, "Requires all arguments to be a floating point type.");

    // For one discussion on these variadic coding idioms, see: https://arne-mertz.de/2016/11/more-variadic-templates/
    (void) std::initializer_list<int>{
        (t < vs ? std::swap(t, vs), 0 : 0)...};

    value_type smallSum(0);
    (void) std::initializer_list<int>{
        (smallSum += std::exp(vs - t),0)... };

    return t + log1p_switch(smallSum);
}
