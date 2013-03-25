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
#include "blt_util/blt_exception.hh"

#ifdef KILL_EXCEPTIONS
#include "blt_util/log.hh"

#include <cstdlib>

#include <iostream>
#endif



blt_exception::
blt_exception(const char* s)
    : message(s) {
#ifdef KILL_EXCEPTIONS
    log_os << "ERROR:: " << s << std::endl;
    abort();
#endif
}

