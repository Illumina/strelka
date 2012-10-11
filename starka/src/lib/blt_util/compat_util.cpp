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

#include "compat_util.hh"

#include <cmath>
#include <cstring>

#include <iostream>


double
compat_round(const double x) {
    if(x>=0.) {
        return std::floor(x+0.5);
    } else {
        return std::ceil(x-0.5);
    }
}



const char*
compat_basename(const char* str) {
#ifdef _WIN32
    static const char pathsep('\\');
#else
    static const char pathsep('/');
#endif
    const char* res(strrchr(str,pathsep));
    if(NULL==res) return str;
    return res+1;
}



static
void
test_round(const double x) {
    std::cerr << "test_round input: " << x << " output: " << compat_round(x) <<"\n";
}

static
void
test_bn(const char* x) {
    std::cerr << "test_bn input: '" << x << "' output: '" << compat_basename(x) <<"'\n";
}

void
test_compat() {
  test_round(3.5);
  test_round(3.2);
  test_round(3.7);
  test_round(-1);
  test_round(-1.5);
  test_round(-1.7);
  test_round(-1.2);

  test_bn("/foo");
  test_bn("foo");
  test_bn("//");
  test_bn("/");
  test_bn("");
}


//int main() { test_compat(); }
