//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

#include "strelka_info.hh"
#include "strelka_option_parser.hh"
#include "strelka_shared.hh"

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
       "usage: " << name() << " [options]\n\n" << visible << "\n";

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
