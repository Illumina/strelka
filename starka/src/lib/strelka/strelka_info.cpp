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

#include "starling_common/starling_option_parser.hh"
#include "strelka/strelka_info.hh"
#include "strelka/strelka_option_parser.hh"
#include "strelka/strelka_shared.hh"

#include "blt_util/log.hh"

#include <cstdlib>

#include <iostream>



void
strelka_info::
usage(const char* xmessage) const {

    std::ostream& os(log_os);

    static strelka_options default_opt;
    static const po::options_description visible(get_strelka_option_parser(default_opt));

    os <<
        "\n" << name() << " - joint tumor/normal small-variant caller\n"
        "\tversion: " << version() << "\n"
        "\n"
        "usage: " << name() << " [options]\n\n" << visible
             << "\n\n[ ***** new single-sample options ***** ]\n\n";

    static const po::options_description visible2(get_starling_shared_option_parser(default_opt));
    os << visible2
       << "\n\n\n[ ***** legacy single-sample options ***** ]\n\n";
    write_starling_legacy_options(os);
    os << "\n";

    if(xmessage) {
        os << "\n"
           << "******** COMMAND-LINE ERROR:: " << xmessage << " ********\n"
           << "\n";
        exit(EXIT_FAILURE);
    }

    exit(EXIT_SUCCESS);
}



void
strelka_info::
doc() const {

    std::ostream& os(log_os);

    os << "\n" << name() << " documentation page TBD.\n";

    exit(EXIT_SUCCESS);
}
