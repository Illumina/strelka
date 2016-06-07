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

/// \author Chris Saunders
///

#pragma once

#include "parse_util.hh"

#include "boost/program_options.hpp"

#include <string>
#include <vector>


/// holds floating-point values with exact control over pretty-print based on inital
/// string input
///
template <typename FLOAT_TYPE>
struct PrettyFloat {
    PrettyFloat(
        const char* init = nullptr)
        : _numval(0)
    {
        if (nullptr == init) {
            _strval="0";
            return;
        }
        update(init);
    }

    void
    update(const char* str) {
        _strval=str;
        _numval=illumina::blt_util::parse_double_str(_strval);
    }

    const std::string&
    strval() const { return _strval; }

    FLOAT_TYPE
    numval() const { return _numval; }

private:
    std::string _strval;
    FLOAT_TYPE _numval;
};


/// this allows PrettyFloatPrinter to work with boost program_options
template <typename FLOAT_TYPE>
void validate(boost::any& v,
              const std::vector<std::string>& values,
              PrettyFloat<FLOAT_TYPE>*, int)
{
    if (v.empty()) {
        v = boost::any(PrettyFloat<FLOAT_TYPE>());
    }
    PrettyFloat<FLOAT_TYPE>* tv = boost::any_cast<PrettyFloat<FLOAT_TYPE>>(&v);
    assert(nullptr != tv);

    if (values.size() != 1) {
        throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
    }

    tv->update(values[0].c_str());
}

