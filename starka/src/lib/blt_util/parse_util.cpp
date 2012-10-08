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

#include "blt_util/parse_util.hh"
#include "common/Exceptions.hh"

#include <boost/format.hpp>

#include <cerrno>
#include <climits>
#include <cstdlib>

#include <limits>


namespace casava
{
namespace blt_util
{


static
void
parse_exception(const char* type_label,
                const char* parse_str) {
    BOOST_THROW_EXCEPTION(
        casava::common::CasavaException(errno,(boost::format("ERROR: Can't parse %s from string: '%s'\n ") % type_label % parse_str).str()));
}



unsigned
parse_unsigned(const char*& s){

    static const int base(10);

    errno = 0;

    char* endptr;
    const unsigned long val(strtoul(s, &endptr, base));
    if ((errno == ERANGE && (val == ULONG_MAX || val == 0))
        || (errno != 0 && val == 0) || (endptr == s)) {
        parse_exception("unsigned long",s);
    }

    if(val > std::numeric_limits<unsigned>::max()){
        parse_exception("unsigned",s);
    }

    s = endptr;

    return static_cast<unsigned>(val);
}



unsigned
parse_unsigned_str(const std::string& s) {
    const char* s2(s.c_str());
    const unsigned val(parse_unsigned(s2));
    if((s2-s.c_str())!=static_cast<int>(s.length())) {
        parse_exception("unsigned",s.c_str());
    }
    return val;
}



int
parse_int(const char*& s){

    static const int base(10);

    errno = 0;

    char* endptr;
    const long val(strtol(s, &endptr, base));
    if ((errno == ERANGE && (val == LONG_MIN || val == LONG_MAX))
        || (errno != 0 && val == 0) || (endptr == s)) {
        parse_exception("long int",s);
    }

    if(val > std::numeric_limits<int>::max()){
        parse_exception("int",s);
    }

    s = endptr;

    return static_cast<int>(val);
}



int
parse_int_str(const std::string& s) {
    const char* s2(s.c_str());
    const int val(parse_int(s2));
    if((s2-s.c_str())!=static_cast<int>(s.length())) {
        parse_exception("int",s.c_str());
    }
    return val;
}



}
}
