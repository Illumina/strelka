// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
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
///
/// \author Chris Saunders
///

#pragma once


#include "starling_common/indel_key.hh"
#include "starling_common/indel_data.hh"

#include <iosfwd>


// Represents a single indel observation
struct indel_observation
{
    indel_key key;
    indel_observation_data data;
};



// Represents all information about an indel
struct indel
{
    indel_key key;
    indel_data data;
};



// Debugging dump:
std::ostream& operator<<(std::ostream& os, const indel_observation& obs);
std::ostream& operator<<(std::ostream& os, const indel& in);

