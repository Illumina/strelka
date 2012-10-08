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

/// \author Chris Saunders
///
#include "blt_util/stream_stat.hh"

#include <iostream>



std::ostream&
operator<<(std::ostream& os,const stream_stat& ss) {

    os << "sample_size: " << ss.sample_size() << " min: " << ss.min() << " max: " << ss.max()
       << " mean: " << ss.mean() << " sd: " << ss.sd() << " se: " << ss.stderror();


    return os;
}
