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

#ifndef __STRELKA_INFO_HH
#define __STRELKA_INFO_HH

#include "blt_util/prog_info.hh"


struct strelka_info : public prog_info {

    static
    const prog_info& get() {
        static const strelka_info vci;
        return vci;
    }

private:
    const char* name() const {
        static const char NAME[] = "strelka";
        return NAME;
    }

    const char* version() const {
        static const char VERSION[] = "${VERSION}";
        return VERSION;
    }

    void usage(const char* xmessage = 0) const;

    void doc() const;

    strelka_info() {}
    virtual ~strelka_info() {}
};

#endif
