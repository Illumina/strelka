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

#ifndef __STRINGER_HH
#define __STRINGER_HH

#include "compat_util.hh"
#include "scan_string.hh"

#include <cstdio>

#include <typeinfo>



struct stringer_base {

    stringer_base()
        : _scanstr(NULL)
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
struct stringer : public stringer_base {

    stringer<T>()
    {
        _scanstr=scan_string<T>();
        if (NULL==_scanstr) { type_error(typeid(T).name()); }
    }


    const char*
    get32(const T val) const {
        static const unsigned buff_size(32);
        const int write_size(snprintf(_buff32,buff_size,_scanstr,val));
        if ((write_size<0) || (write_size >= static_cast<int>(buff_size))) {
            get32_error(write_size);
        }
        return _buff32;
    }
};


#endif
