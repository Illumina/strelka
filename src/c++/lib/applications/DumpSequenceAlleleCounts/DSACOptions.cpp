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

#include "DSACOptions.hh"

#include "blt_util/log.hh"
#include "common/ProgramUtil.hh"

#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"



static
void
usage(
    std::ostream& os,
    const illumina::Program& prog,
    const boost::program_options::options_description& visible,
    const char* msg = nullptr)
{
    usage(os, prog, visible, "summarize binary allele counts file to stdout", " > counts_dump", msg);
}



void
parseDSACOptions(
    const illumina::Program& prog,
    int argc,
    char** argv,
    DSACOptions& opt)
{
    namespace po = boost::program_options;
    po::options_description req("configuration");
    req.add_options()
    ("exclude-basecalls", po::value(&opt.isExcludeBasecalls)->zero_tokens(),
     "Exclude basecall info from the dump")
    ("exclude-indels", po::value(&opt.isExcludeIndels)->zero_tokens(),
     "Exclude indel info from the dump")
    ("counts-file", po::value(&opt.countsFilename),
     "read binary allele counts from filename (required, no default)")
    ("extended", po::value(&opt.isExtendedOutput)->zero_tokens(),
     "print extended output in tab-delimited format, which should be more suitable for analysis with external tools")
    ;

    po::options_description help("help");
    help.add_options()
    ("help,h","print this message");

    po::options_description visible("options");
    visible.add(req).add(help);

    bool po_parse_fail(false);
    po::variables_map vm;
    try
    {
        po::store(po::parse_command_line(argc, argv, visible,
                                         po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
        po::notify(vm);
    }
    catch (const boost::program_options::error& e)     // todo:: find out what is the more specific exception class thrown by program options
    {
        log_os << "\nERROR: Exception thrown by option parser: " << e.what() << "\n";
        po_parse_fail=true;
    }

    if ((argc<=1) || (vm.count("help")) || po_parse_fail)
    {
        usage(log_os,prog,visible);
    }

    if (opt.countsFilename.empty())
    {
        usage(log_os,prog,visible,"Must specify counts file");
    }
    if (! boost::filesystem::exists(opt.countsFilename))
    {
        usage(log_os,prog,visible,"Counts file does not exist");
    }
}
