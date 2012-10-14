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

#include "blt_exception.hh"
#include "parse_util.hh"

#include <cerrno>
#include <climits>
#include <cstdlib>

#include <limits>
#include <sstream>



static
void
parse_exception(const char* type_label,
                const char* parse_str) {

    std::ostringstream oss;
    oss << "ERROR: Can't parse " << type_label << " from string: '" << parse_str << "'";
    throw blt_exception(oss.str().c_str());
}


namespace casava {
namespace blt_util {

unsigned
parse_unsigned(const char*& s){

    static const int base(10);

    errno = 0;

    char* endptr;
    const unsigned long val(strtoul(s, &endptr, base));
    if ((errno == ERANGE && (val == ULONG_MAX || val == 0))
        || (errno != 0 && val == 0) || (endptr == s)) {
        parse_exception("unsigned long",s);
    }

    if(val > std::numeric_limits<unsigned>::max()){
        parse_exception("unsigned",s);
    }

    s = endptr;

    return static_cast<unsigned>(val);
}



unsigned
parse_unsigned_str(const std::string& s) {
    const char* s2(s.c_str());
    const unsigned val(parse_unsigned(s2));
    if((s2-s.c_str())!=static_cast<int>(s.length())) {
        parse_exception("unsigned",s.c_str());
    }
    return val;
}



int
parse_int(const char*& s){

    const char* endptr(s);
    const long val(parse_long(endptr));

    if(val > std::numeric_limits<int>::max()){
        parse_exception("int",s);
    }

    s = endptr;

    return static_cast<int>(val);
}



int
parse_int_str(const std::string& s) {
    const char* s2(s.c_str());
    const int val(parse_int(s2));
    if((s2-s.c_str())!=static_cast<int>(s.length())) {
        parse_exception("int",s.c_str());
    }
    return val;
}



long
parse_long(const char*& s){

    static const int base(10);

    errno = 0;

    char* endptr;
    const long val(strtol(s, &endptr, base));
    if ((errno == ERANGE && (val == LONG_MIN || val == LONG_MAX))
        || (errno != 0 && val == 0) || (endptr == s)) {
        parse_exception("long int",s);
    }

    s = endptr;

    return val;
}



long
parse_long_str(const std::string& s) {
    const char* s2(s.c_str());
    const long val(parse_long(s2));
    if((s2-s.c_str())!=static_cast<int>(s.length())) {
        parse_exception("long int",s.c_str());
    }
    return val;
}



double
parse_double(const char*& s){

    errno = 0;

    char* endptr;
    const double val(strtod(s, &endptr));
    if ((errno == ERANGE) || (endptr == s)) {
        parse_exception("double",s);
    }

    s = endptr;
    return val;
}



double
parse_double_str(const std::string& s) {
    const char* s2(s.c_str());
    const double val(parse_double(s2));
    if((s2-s.c_str())!=static_cast<int>(s.length())) {
        parse_exception("double",s.c_str());
    }
    return val;
}
 
}
}
