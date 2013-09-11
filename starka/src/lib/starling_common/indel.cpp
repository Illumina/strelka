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
///
/// \author Chris Saunders
///

#include "starling_common/indel.hh"

#include <iostream>



std::ostream&
operator<<(std::ostream& os,
           const indel_observation& obs) {

    os << obs.key << obs.data;

    return os;
}

std::ostream&
operator<<(std::ostream& os,
           const indel& in) {

    os << in.key << in.data;

    return os;
}
