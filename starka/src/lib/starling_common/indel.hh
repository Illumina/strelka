// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
//

/// \file
///
/// \author Chris Saunders
///

#ifndef __INDEL_HH
#define __INDEL_HH

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
std::ostream& operator<<(std::ostream& os, const indel& in);



#endif
