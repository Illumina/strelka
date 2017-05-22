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

#include "GetSequenceErrorCounts.hh"
#include "GetSequenceErrorCountsInfo.hh"
#include "GetSequenceErrorCountsRun.hh"
#include "SequenceErrorCountsOptionsParser.hh"
#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"
#include "common/Exceptions.hh"
#include "starling_common/starling_arg_parse.hh"

#include <cassert>
#include <cstdlib>



namespace
{
const prog_info& pinfo(GetSequenceErrorCountsInfo::get());
}



void
GetSequenceErrorCounts::
runInternal(int argc, char* argv[]) const
{
    SequenceErrorCountsOptions opt;

    for (int i(0); i<argc; ++i)
    {
        if (i) opt.cmdline += ' ';
        opt.cmdline += argv[i];
    }

    std::vector<std::string> legacy_starling_args;
    po::variables_map vm;
    try
    {
        po::options_description visible(getSequenceErrorCountsOptionsParser(opt));
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

    finalizeSequenceErrorCountsOptions(pinfo,vm,opt);

    getSequenceErrorCountsRun(pinfo,opt);
}
