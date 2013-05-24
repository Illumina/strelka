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

#pragma once


#include "starling_common/indel_key.hh"
#include "starling_common/indel_data.hh"

#include <iosfwd>


// Represents a single indel observation
struct indel_observation {
    indel_key key;
    indel_observation_data data;
};



// Represents all information about an indel
struct indel {
    indel_key key;
    indel_data data;
};



// Debugging dump:
std::ostream& operator<<(std::ostream& os, const indel_observation& obs);
std::ostream& operator<<(std::ostream& os, const indel& in);

