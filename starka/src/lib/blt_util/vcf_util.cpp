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


#include "blt_util/vcf_util.hh"

#include <cassert>
#include <ctime>

#include <iostream>



std::ostream&
vcf_fileDate(std::ostream& os) {
    const time_t t(time(NULL));
    struct tm* ct(localtime(&t));
    assert(NULL != ct);

    static const unsigned dsize(64);
    char datebuf[dsize];
    const size_t ret(strftime(datebuf,dsize,"%Y%m%d",ct));
    assert(ret!=0);
    return os << datebuf;
}
