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
///
/// \author Chris Saunders
///

#include "starling_common/indel_key.hh"

#include <iostream>



std::ostream&
operator<<(std::ostream& os,
           const indel_key& ik) {
    os << "INDEL pos: " << ik.pos
       << " type: " << INDEL::get_index_label(ik.type)
       << " len: " << ik.length
       << " swap_dlen: " << ik.swap_dlength << "\n";
    return os;
}


