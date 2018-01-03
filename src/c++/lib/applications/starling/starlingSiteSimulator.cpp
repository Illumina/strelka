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

///
/// \author Chris Saunders
///
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the beginning of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///


#include "starlingSiteSimulator.hh"

#include "starling_sim_test.hh"

#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"
#include "common/Exceptions.hh"

#include "boost/program_options.hpp"



void
starlingSiteSimulator::
runInternal(int argc,char* argv[]) const
{
    starling_options opt;
    starling_site_sim_options sim_opt;

    for (int i(0); i<argc; ++i)
    {
        if (i) opt.cmdline += ' ';
        opt.cmdline += argv[i];
    }

    namespace po = boost::program_options;
    po::options_description req("configuration");
    req.add_options()
    ("total-sites", po::value(&sim_opt.total_sites)->default_value(sim_opt.total_sites),"number of sites to simulate")
    ("cov", po::value(&sim_opt.coverage)->default_value(sim_opt.coverage),"site sequencing depth")
    ("exact-cov","By default simulator uses Poisson coverage distribution at requested mean, this flag changes to exact coverage")
    ("variants-only","only simulate non-reference sites")
    ("ref-only","only simulate reference sites")
    ("het-only","all variants are hets")
    ("seed",po::value<uint32_t>(&sim_opt.seed),"seed")
    ("qscores",po::value(&sim_opt.qval_file),"tab-delimited file specifying basecall qscore distribution (default: all basecalls are Q30)")
    ;

    po::options_description visible("options");
    visible.add(req);

    bool po_parse_fail(false);
    po::variables_map vm;
    try
    {
        po::store(po::parse_command_line(argc, argv, visible), vm);
        po::notify(vm);
    }
    catch (const boost::program_options::error& e)
    {
        log_os << "\nERROR: Exception thrown by option parser: " << e.what() << "\n";
        po_parse_fail=true;
    }

    if ((argc<=1) || (vm.count("help")) || po_parse_fail)
    {
        log_os << "\n starling site simulator...\n\n";
        log_os << "usage: program [options] > called\n\n";
        log_os << visible << "\n";
        exit(EXIT_FAILURE);
    }

    if (vm.count("variants-only"))
    {
        sim_opt.mode=SIM_GERMLINE;
    }

    if (vm.count("ref-only"))
    {
        sim_opt.mode=SIM_REF;
    }

    if (vm.count("het-only"))
    {
        sim_opt.is_het_only = true;
    }

    if (vm.count("exact-cov"))
    {
        sim_opt.is_exact_cov = true;
    }

    starling_site_sim(opt,sim_opt);
}
