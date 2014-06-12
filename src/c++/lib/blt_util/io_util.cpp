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

/// \author Chris Saunders
///

#include "blt_util/io_util.hh"
#include "blt_util/log.hh"

#include <cstdlib>

#include <fstream>
#include <iostream>



void
open_ifstream(std::ifstream& ifs,
              const char* filename) {

    ifs.open(filename);
    if (! ifs) {
        log_os << "ERROR: Can't open file: " << filename << "\n";
        exit(EXIT_FAILURE);
    }
}
