// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

/// \file

/// \author Chris Saunders
///

#include "compat_util.hh"

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <iostream>


#ifdef _WIN32
#include "win32_realpath.c"
#endif



bool
compat_realpath(std::string& path) {
    errno=0;
    const char* newpath(realpath(path.c_str(),NULL));
    if((NULL==newpath) || (errno!=0)) {
        if(NULL!=newpath) free((void*)newpath);
        return false;
    }
    path = newpath;
    free((void*)newpath);
    return true;
}



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
