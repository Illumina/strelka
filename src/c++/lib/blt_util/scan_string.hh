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

/// \author Chris Saunders
///

#ifndef __SCAN_STRING_HH
#define __SCAN_STRING_HH

#include "boost/static_assert.hpp"


template <typename T>
const char*
scan_string()
{
    // no scan_string available for type:
    BOOST_STATIC_ASSERT(sizeof(T)==0);
    return NULL;
}


template <>
inline
const char*
scan_string<int>()
{
    return "%d";
}

template <>
inline
const char*
scan_string<long>()
{
    return "%ld";
}


#endif
