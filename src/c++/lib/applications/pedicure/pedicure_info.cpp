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

///
/// \author Chris Saunders
///

#include "inovo_info.hh"
#include "inovo_option_parser.hh"
#include "inovo_shared.hh"

#include "blt_util/log.hh"
#include "starling_common/starling_base_option_parser.hh"

#include <cstdlib>

#include <iostream>



void
inovo_info::
usage(const char* xmessage) const
{
    std::ostream& os(log_os);

    static inovo_options default_opt;
    static const po::options_description visible(get_inovo_option_parser(default_opt));

    os <<
       "\n" << name() << " - de-novo small-variant caller\n"
       "\tversion: " << version() << "\n"
       "\n"
       "usage: " << name() << " [options]\n\n" << visible
       << "\n\n\n[ ***** legacy single-sample options ***** ]\n\n";
    write_starling_legacy_options(default_opt,os);
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
inovo_info::
doc() const
{

    std::ostream& os(log_os);

    os << "\n" << name() << " documentation page TBD.\n";

    exit(EXIT_SUCCESS);
}
