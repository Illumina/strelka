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

#include "Tier2OptionsParser.hh"
#include "blt_common/blt_arg_validate.hh"



boost::program_options::options_description
getTier2OptionsDescription(
    Tier2Options& opt)
{
    namespace po = boost::program_options;
    po::options_description desc("Tier2 data thresholds");
    desc.add_options()
    ("tier2-min-mapping-quality",
     po::value(&opt.minMappingErrorPhredProb),
     "Min mapping quality used for tier2 calling")
    ;
    return desc;
}



bool
parseTier2Options(
    const boost::program_options::variables_map& /*vm*/,
    Tier2Options& /*opt*/,
    std::string& errorMsg)
{
    errorMsg.clear();

    // no extra parsing steps apart from what boost::program_options is already doing

    return (! errorMsg.empty());
}
