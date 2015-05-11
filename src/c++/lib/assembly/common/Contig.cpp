// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Rumovsky
// Copyright (c) 2013-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Ole Schulz-Trieglaff
///

#include "Contig.hh"

#include <iostream>


std::ostream&
operator<<(std::ostream& os, const Contig& contig)
{
    os << ">CTG::LEN=" << contig.seq.size()
       << ":SEEDS="    << contig.seedReadCount
       << ":COV="      << contig.avgCoverage << "\n"
       << contig.seq;
    //os << "\n";
    return os;
}
