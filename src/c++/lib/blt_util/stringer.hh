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

/// \author Chris Saunders
///

#pragma once

#include "compat_util.hh"
#include "scan_string.hh"

#include <cstdio>

#include <typeinfo>


struct stringer_base
{
    stringer_base()
        : _scanstr(nullptr)
    {}

protected:
    static
    void type_error(const char* tiname);
    void get32_error(const int write_size) const;

    mutable char _buff32[32];
    const char* _scanstr;
};


/// String conversion utility which is harder-to-use but faster than stringstream/lexical_cast
///
/// Safety notes:
/// 1) client must create one object for each thread
/// 2) The string pointer returned will be invalid at the next conversion call to stringer
///
template <typename T>
struct stringer : public stringer_base
{
    stringer<T>()
    {
        _scanstr=scan_string<T>();
        if (nullptr==_scanstr)
        {
            type_error(typeid(T).name());
        }
    }

    const char*
    get32(const T val) const
    {
        static const unsigned buff_size(32);
        const int write_size(snprintf(_buff32,buff_size,_scanstr,val));
        if ((write_size<0) || (write_size >= static_cast<int>(buff_size)))
        {
            get32_error(write_size);
        }
        return _buff32;
    }
};
