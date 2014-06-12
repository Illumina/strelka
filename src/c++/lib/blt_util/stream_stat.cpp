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

/// \author Chris Saunders
///
#include "stream_stat.hh"

#include <iostream>



std::ostream&
operator<<(std::ostream& os,const stream_stat& ss) {

    os << "sample_size: " << ss.size() << " min: " << ss.min() << " max: " << ss.max()
       << " mean: " << ss.mean() << " sd: " << ss.sd() << " se: " << ss.stderror();


    return os;
}
