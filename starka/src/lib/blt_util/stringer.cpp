// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
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
