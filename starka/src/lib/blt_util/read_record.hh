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
#ifndef __READ_RECORD_HH
#define __READ_RECORD_HH


#include "boost/utility.hpp"


struct read_record : private boost::noncopyable {

    virtual ~read_record() {}

    virtual const char* qname() const = 0;

    virtual bool is_unmapped() const = 0;
    virtual bool is_dup() const = 0;
    virtual bool is_filter() const = 0;

    virtual int pos() const = 0;
};


#endif
