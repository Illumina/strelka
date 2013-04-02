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
