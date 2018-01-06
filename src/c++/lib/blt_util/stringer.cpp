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

#include "blt_exception.hh"
#include "stringer.hh"

#include <sstream>



void
stringer_base::
type_error(const char* tiname)
{
    std::ostringstream oss;
    oss << "Can't initialize stringer object for type: "  << tiname;
    throw blt_exception(oss.str().c_str());
}



void
stringer_base::
get32_error(const int write_size) const
{
    _buff32[31]='\0';
    std::ostringstream oss;
    oss << "stringer.get32 failed to convert type to string. write_size: '" << write_size << " buffer: " << _buff32;
    throw blt_exception(oss.str().c_str());
}
