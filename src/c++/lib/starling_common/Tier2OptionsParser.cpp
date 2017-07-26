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
     po::value(&opt.tier2_min_mapping_quality),
     "Activate tier2 call validation and use the following mapping quality threshold in the tier2 set")
    ("tier2-mismatch-density-filter-count",
     po::value(&opt.tier2_mismatch_density_filter_count),
     "Use the specified less stringent mismatch density count in tier2 evaluation.")
    ("tier2-no-mismatch-density-filter",
     po::value(&opt.is_tier2_no_mismatch_density_filter)->zero_tokens(),
     "Don't apply mismatch density filter to tier2 data.")
    ("tier2-include-singleton",
     po::value(&opt.is_tier2_include_singleton)->zero_tokens(),
     "Don't filter singleton read pairs from tier2 data (will have no effect without rescue mode).")
    ("tier2-include-anomalous",
     po::value(&opt.is_tier2_include_anomalous)->zero_tokens(),
     "Don't filter anomalous read pairs from tier2 data (will probably have no effect without rescue mode).")
    ("tier2-indel-nonsite-match-prob",
     po::value(&opt.tier2_indel_nonsite_match_prob),
     "Reset indel-nonsite-match-prob for tier2 analysis.")
    ;
    return desc;
}



bool
parseTier2Options(
    const boost::program_options::variables_map& vm,
    Tier2Options& opt,
    std::string& errorMsg)
{
    errorMsg.clear();

    if (vm.count("tier2-indel-nonsite-match-prob"))
    {
        opt.is_tier2_indel_nonsite_match_prob=true;
    }

    if (vm.count("tier2-min-mapping-quality"))
    {
        opt.is_tier2_min_mapping_quality=true;
    }

    if (vm.count("tier2-mismatch-density-filter-count"))
    {
        opt.is_tier2_mismatch_density_filter_count=true;
    }

    check_option_arg_range(opt.tier2_indel_nonsite_match_prob,"tier2-indel-nonsite-match-prob",0.,1.,errorMsg);
    return (! errorMsg.empty());
}
