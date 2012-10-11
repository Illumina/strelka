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
#include "stringer.hh"

#include <sstream>



void
stringer_base::
type_error(const char* tiname) {
    std::ostringstream oss;
    oss << "ERROR: Can't initialize stringer object for type: "  << tiname << "\n";
    throw blt_exception(oss.str().c_str());
}



void
stringer_base::
get32_error(const int write_size) const {
    _buff32[31]='\0';
    std::ostringstream oss;
    oss << "ERROR: stringer.get32 failed to convert type to string. write_size: '" << write_size << " buffer: " << _buff32 << "'\n";
    throw blt_exception(oss.str().c_str());
}
