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

#include "starling_common/pos_basecall_buffer.hh"

#include <iostream>



// debug dumpers:
void
pos_basecall_buffer::
dump(std::ostream& os) const {

    pciter i(_pdata.begin()), i_end(_pdata.end());
    for(; i!=i_end; ++i) {
        os << "pc_buff pos: " << i->first << "\n";
    }
}
