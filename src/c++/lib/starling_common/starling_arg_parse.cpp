//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
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

#include "blt_common/blt_arg_parse_util.hh"
#include "starling_common/starling_arg_parse.hh"

//#define DEBUG_PARSER

#ifdef DEBUG_PARSER
#include "blt_util/log.hh"
#endif


void
legacy_starling_arg_parse(
    arg_data& ad,
    starling_base_options& opt)
{
    const prog_info& pinfo(ad.pinfo);

    if ((ad.size()==1) ||
        ((ad.size()==2) && ((ad.argstr[1] == "-h") || (ad.argstr[1] == "-help") ||
                            (ad.argstr[1] == "--help") || (ad.argstr[1] == "-")))) pinfo.usage();

    opt.cmdline = ad.cmdline;

    bool is_max_indel_size(false);

    const unsigned as(ad.size());
    for (unsigned i(0); i<as; ++i)
    {
        if (ad.argmark[i]) continue;

        if (ad.argstr[i]=="-max-indel-size")
        {
            set_arg(i,ad,is_max_indel_size,opt.maxIndelSize);
        }
        else if (ad.argstr[i]=="-h")
        {
            pinfo.usage();
        }
        else
        {
            continue;
        }

        ad.argmark[i] = true;
    }

    ad.finalize_args();
}
