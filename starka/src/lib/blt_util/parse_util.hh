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

#ifndef __PARSE_UTIL_HH
#define __PARSE_UTIL_HH

#include <ciso646>

#include <string>


namespace casava
{
namespace blt_util
{

/// parse an unsigned from char* with several error checks, and advance
/// pointer to end of unsigned input
///
unsigned
parse_unsigned(const char*& s);

/// std::string version, no ptr advance obviously. explicit rename
/// gaurds against unexpected std::string temporaries
unsigned
parse_unsigned_str(const std::string& s);


/// parse an int from char* with several error checks, and advance
/// pointer to end of integer input
///
int
parse_int(const char*& s);

/// std::string version, no ptr advance obviously. explicit rename
/// gaurds against unexpected std::string temporaries
int
parse_int_str(const std::string& s);

}
}

#endif
