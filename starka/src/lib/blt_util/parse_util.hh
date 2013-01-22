// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use, copy,
// modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//

/// \file

/// \author Chris Saunders
///

#ifndef __PARSE_UTIL_HH
#define __PARSE_UTIL_HH

#include "boost/static_assert.hpp"

#include <string>


namespace casava {
namespace blt_util {

/// parse TYPE from char* with several error checks, and advance
/// pointer to end of TYPE input
///
unsigned
parse_unsigned(const char*& s);

int
parse_int(const char*& s);

long
parse_long(const char*& s);

double
parse_double(const char*& s);



/// std::string version of above, no ptr advance obviously. explicit rename
/// of functions guards against unexpected std::string temporaries
///
unsigned
parse_unsigned_str(const std::string& s);

int
parse_int_str(const std::string& s);

long
parse_long_str(const std::string& s);

double
parse_double_str(const std::string& s);



/// template version:
///
template <typename T>
T
parse_type(const char*&) {
    // no scan_string available for type:
    BOOST_STATIC_ASSERT(sizeof(T)==0);
    return T();
}


template <>
inline
unsigned
parse_type<unsigned>(const char*& s) { return parse_unsigned(s); }

template <>
inline
int
parse_type<int>(const char*& s) { return parse_int(s); }

template <>
inline
long
parse_type<long>(const char*& s) { return parse_long(s); }

template <>
inline
double
parse_type<double>(const char*& s) { return parse_double(s); }


}
}

#endif
