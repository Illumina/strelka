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

#ifndef __PROG_INFO_HH
#define __PROG_INFO_HH

struct prog_info {
    virtual
    const char* name() const = 0;

    virtual
    const char* version() const = 0;

    virtual
    void usage(const char* xmessage = 0) const = 0;

    virtual
    void doc() const = 0;

protected:
    prog_info() {}
    virtual ~prog_info() {}
};

#endif
