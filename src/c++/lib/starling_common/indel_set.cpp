// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

/// \file
///
/// \author Chris Saunders
///

#include "starling_common/indel_set.hh"

#include <iostream>



void
dump_indel_set(const indel_set_t& is,
               std::ostream& os) {

    indel_set_t::const_iterator i(is.begin()), i_end(is.end());
    for (; i!=i_end; ++i) os << *i;
}
