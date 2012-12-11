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

#include "blt_util/log.hh"
#include "blt_util/sig_handler.hh"

#include <cstdlib>
#include <signal.h>

#include <iostream>
#include <string>


std::string _progname;
std::string _cmdline;



void
blt_sig_handler (int sig) {

    switch (sig) {
    case SIGTERM:
        log_os << "ERROR: " << _progname << " received termination signal. cmdline: " << _cmdline << std::endl;
        exit(EXIT_FAILURE);
#ifndef _WIN32
    case SIGINT:
        log_os << "ERROR: " << _progname << " received interupt signal. cmdline: " << _cmdline << std::endl;
        exit(EXIT_FAILURE);
#endif
    default:
        log_os << "INFO: " << _progname << " received signal no: " << sig << std::endl;
        break;
    }
}



void
initialize_blt_signals(const char* progname,
                       const char* cmdline) {

    _progname=progname;
    _cmdline=cmdline;

    signal(SIGTERM, blt_sig_handler);
#ifndef _WIN32
    signal(SIGINT, blt_sig_handler);
#endif
}
