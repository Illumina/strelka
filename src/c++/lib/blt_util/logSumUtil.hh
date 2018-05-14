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
/// \brief Various utilities for taking sums in log-space
///

#pragma once

#include "blt_util/math_util.hh"

#include "boost/concept_check.hpp"


/// Returns the equivalent of log(exp(x1)+exp(x2))
///
template <typename FloatType>
FloatType
getLogSum(FloatType x1, FloatType x2)
{
    static_assert(std::is_floating_point<FloatType>::value, "Requires floating point type.");

    if (x1<x2) std::swap(x1,x2);
    return x1 + log1p_switch(std::exp(x2-x1));
}

/// Returns the equivalent of log(exp(x1)+exp(x2)+exp(x3))
///
template <typename FloatType>
FloatType
getLogSum(FloatType x1, FloatType x2, FloatType x3)
{
    static_assert(std::is_floating_point<FloatType>::value, "Requires floating point type.");

    if (x1<x2) std::swap(x1,x2);
    if (x1<x3) std::swap(x1,x3);
    return x1 + log1p_switch(std::exp(x2-x1)+std::exp(x3-x1));
}

/// Returns the equivalent of log(exp(x1)+exp(x2)+exp(x3)+exp(x4))
///
template <typename FloatType>
FloatType
getLogSum(FloatType x1, FloatType x2, FloatType x3, FloatType x4)
{
    static_assert(std::is_floating_point<FloatType>::value, "Requires floating point type.");

    if (x1<x2) std::swap(x1,x2);
    if (x1<x3) std::swap(x1,x3);
    if (x1<x4) std::swap(x1,x4);
    return x1 + log1p_switch(std::exp(x2-x1)+std::exp(x3-x1)+std::exp(x4-x1));
}


/// Returns the equivalent of log(exp(x_1)+exp(x_2)+....exp(x_n))
///
/// This is the iterator version of getLogSum, it is slightly less fast (~20%) than the version coded for a specific
/// number of arguments, even before considering any runtime penalties to instantiating the iterated container.
///
template <typename IterType>
typename std::iterator_traits<IterType>::value_type
getLogSumSequence(
    IterType beginIter,
    IterType endIter)
{
    // Concept assertions removed due to unused typedef warning in clang & boost 1.58, this is supposed to be fixed in boost 1.59
    //
    //BOOST_CONCEPT_ASSERT((boost::ForwardIterator<IterType>));
    typedef typename std::iterator_traits<IterType>::value_type value_type;
    static_assert(std::is_floating_point<value_type>::value, "Requires iterator on floating point type.");

    static const value_type neginf(-std::numeric_limits<value_type>::infinity());

    if (beginIter == endIter) return neginf;
    if (beginIter+1 == endIter) return *beginIter;

    IterType largest(beginIter);
    for (IterType iter(beginIter+1); iter != endIter; ++iter)
    {
        if (*iter > *largest) largest = iter;
    }

    value_type smallSum(0);
    for (; beginIter != endIter; ++beginIter)
    {
        if (beginIter == largest) continue;
        smallSum += std::exp(*beginIter - *largest);
    }

    return *largest + log1p_switch(smallSum);
}


/// Returns the equivalent of log(exp(x_1)+exp(x_2)+....exp(x_n)), where x_1..x_n are all elements in \p container
///
template <typename ContainerType>
auto
getLogSumSequence(
    const ContainerType& container) -> typename std::iterator_traits<decltype(std::begin(container))>::value_type
{
    return getLogSumSequence(std::begin(container), std::end(container));
}
