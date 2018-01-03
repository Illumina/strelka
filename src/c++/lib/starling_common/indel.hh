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

///
/// \author Chris Saunders
///

#pragma once


#include "starling_common/IndelKey.hh"
#include "starling_common/IndelData.hh"

#include <iosfwd>


/// Represents a single indel observation
struct IndelObservation
{
    IndelKey key;
    IndelObservationData data;
};



/// Represents all information about an indel
struct indel
{
    IndelKey key;
    IndelData data;
};



// Debugging dump:
std::ostream& operator<<(std::ostream& os, const IndelObservation& obs);
std::ostream& operator<<(std::ostream& os, const indel& in);

