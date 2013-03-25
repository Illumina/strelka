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
