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

#include "starling_common/grouper_contig.hh"

#include <iostream>



std::ostream&
operator<<(std::ostream& os,
           const grouper_contig& ctg) {

    os << "GROUPER_CONTIG id: " << ctg.id
       << " chrom: " << ctg.chrom
       << " usable: " << ctg.is_usable << "\n"
       << static_cast<const alignment&>(ctg)
       << "seq: " << ctg.seq << "\n";

    return os;
}
