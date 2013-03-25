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
#ifndef __BLT_EXCEPTION_HH
#define __BLT_EXCEPTION_HH

#include <exception>
#include <string>

/// \brief a minimal exception class
struct blt_exception : public std::exception {

    blt_exception(const char* s);

    ~blt_exception() throw() {}

    const char* what() const throw() { return message.c_str(); }

    std::string message;
};


#endif

