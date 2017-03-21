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

///
/// \author Chris Saunders
///

#include "strelkaNoiseExtractor.hh"
#include "snoise_info.hh"
#include "snoise_run.hh"
#include "snoise_option_parser.hh"

#include "blt_common/blt_arg_parse_util.hh"
#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"
#include "common/Exceptions.hh"
#include "starling_common/starling_arg_parse.hh"
#include <cassert>
#include <cstdlib>
#include "../../starling_common/starling_base_option_parser.hh"


namespace
{
const prog_info& pinfo(snoise_info::get());
}



void
StrelkaNoiseExtractor::
runInternal(int argc, char* argv[]) const
{
    snoise_options opt;

    // set command-line defaults for starling only:
    opt.bsnp_ssd_no_mismatch = 0.35;
    opt.bsnp_ssd_one_mismatch = 0.6;
    opt.max_win_mismatch = 2;
    opt.max_win_mismatch_flank_size = 20;
    opt.is_min_vexp = true;
    opt.min_vexp = 0.25;

    for (int i(0); i<argc; ++i)
    {
        if (i) opt.cmdline += ' ';
        opt.cmdline += argv[i];
    }

    std::vector<std::string> legacy_starling_args;
    po::variables_map vm;
    try
    {
        po::options_description visible(get_snoise_option_parser(opt));
        po::options_description visible2(get_starling_base_option_parser(opt));
        visible.add(visible2);
        po::parsed_options parsed(po::command_line_parser(argc,argv).options(visible).allow_unregistered().run());
        po::store(parsed,vm);
        po::notify(vm);

        // allow remaining options to be parsed using old starling command-line parser:
        legacy_starling_args = po::collect_unrecognized(parsed.options,po::include_positional);
    }
    catch (const boost::program_options::error& e)
    {
        pinfo.usage(e.what());
    }

    if ((argc==1) || vm.count("help"))
    {
        pinfo.usage();
    }

    // temp workaround for blt/starling options which are not (yet)
    // under program_options control:
    //
    arg_data ad(legacy_starling_args,pinfo,opt.cmdline);
    legacy_starling_arg_parse(ad,opt);

    finalize_starling_base_options(pinfo,vm,opt);

    snoise_run(pinfo,opt);
}

