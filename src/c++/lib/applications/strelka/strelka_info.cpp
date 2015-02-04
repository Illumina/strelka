// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
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

#include "strelka_info.hh"
#include "strelka_option_parser.hh"
#include "strelka_shared.hh"
#include "starling_common/starling_base_option_parser.hh"

#include "blt_util/log.hh"

#include <cstdlib>

#include <iostream>



void
strelka_info::
usage(const char* xmessage) const
{
    std::ostream& os(log_os);

    static strelka_options default_opt;
    static const po::options_description visible(get_strelka_option_parser(default_opt));

    os <<
       "\n" << name() << " - joint tumor/normal small-variant caller\n"
       "\tversion: " << version() << "\n"
       "\n"
       "usage: " << name() << " [options]\n\n" << visible
       << "\n\n\n[ ***** legacy single-sample options ***** ]\n\n";
    write_starling_legacy_options(default_opt, os);
    os << "\n";

    if (xmessage)
    {
        os << "\n"
           << "******** COMMAND-LINE ERROR:: " << xmessage << " ********\n"
           << "\n";
        exit(EXIT_FAILURE);
    }

    exit(EXIT_SUCCESS);
}



void
strelka_info::
doc() const
{

    std::ostream& os(log_os);

    os << "\n" << name() << " documentation page TBD.\n";

    exit(EXIT_SUCCESS);
}
