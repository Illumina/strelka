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

/// \author Chris Saunders
///

#pragma once

#include "parse_util.hh"

#include "boost/program_options.hpp"

#include <string>
#include <vector>


/// Holds floating-point values with exact control over pretty-print based on initial string input
///
/// Designed to prevent writing '0.499999999999999999999' where we entered '0.5', in contexts
/// where the output is being read by humans.
///
template <typename FLOAT_TYPE>
struct PrettyFloat
{
    PrettyFloat(
        const char* init = nullptr)
        : _numval(0)
    {
        if (nullptr == init)
        {
            _strval="0";
            return;
        }
        update(init);
    }

    void
    update(const char* str)
    {
        _strval=str;
        _numval=illumina::blt_util::parse_double_str(_strval);
    }

    const std::string&
    strval() const
    {
        return _strval;
    }

    FLOAT_TYPE
    numval() const
    {
        return _numval;
    }

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
    if (v.empty())
    {
        v = boost::any(PrettyFloat<FLOAT_TYPE>());
    }
    PrettyFloat<FLOAT_TYPE>* tv = boost::any_cast<PrettyFloat<FLOAT_TYPE>>(&v);
    assert(nullptr != tv);

    if (values.size() != 1)
    {
        throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
    }

    tv->update(values[0].c_str());
}
