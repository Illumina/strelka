// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
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

    if (xmessage) {
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
