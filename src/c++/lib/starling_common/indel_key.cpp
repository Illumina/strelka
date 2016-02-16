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

#include "starling_common/indel_key.hh"

#include <iostream>

//void
//indel_key::addRanksumInfo(const int mapq, const int baseq, bool is_alt){
//    this->mapq_val=mapq;
//    this->baseq_val=baseq;
//    is_alt = 1; //TODO use alt info
//}


std::ostream&
operator<<(std::ostream& os,
           const indel_key& ik)
{
    os << "INDEL pos: " << ik.pos
       << " type: " << INDEL::get_index_label(ik.type)
       << " len: " << ik.length
       << " mapq: " << ik.mapq_val
       << " baseq: " << ik.baseq_val
       << " swap_dlen: " << ik.swap_dlength << "\n";
    return os;
}


