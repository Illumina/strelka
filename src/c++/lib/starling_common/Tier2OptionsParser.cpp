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

#include "Tier2OptionsParser.hh"
#include "blt_common/blt_arg_validate.hh"



boost::program_options::options_description
getTier2OptionsDescription(
    Tier2Options& opt)
{
    namespace po = boost::program_options;
    po::options_description desc("Tier2 data thresholds");
    desc.add_options()
    ("tier2-min-single-align-score",
     po::value(&opt.tier2_min_single_align_score),
     "Activate tier2 call validation and use the following single alignment score threshold in the tier2 set")
    ("tier2-min-paired-align-score",
     po::value(&opt.tier2_min_paired_align_score),
     "Activate tier2 call validation and use the following paired alignment score threshold in the tier2 set")
    ("tier2-single-align-score-rescue-mode",
     po::value(&opt.is_tier2_single_align_score_rescue_mode)->zero_tokens(),
     "Include non SE-failed reads in tier2 even when a paired score is present and the read is PE-failed")
    ("tier2-mismatch-density-filter-count",
     po::value(&opt.tier2_mismatch_density_filter_count),
     "Use the specified less stringent mismatch density count in tier2 evaluation.")
    ("tier2-no-mismatch-density-filter",
     po::value(&opt.is_tier2_no_mismatch_density_filter)->zero_tokens(),
     "Don't apply mismatch density filter to tier2 data.")
    ("tier2-no-filter-unanchored",
     po::value(&opt.is_tier2_no_filter_unanchored)->zero_tokens(),
     "Don't apply unanchored pair filtration to tier2 data.")
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

    if (vm.count("tier2-min-single-align-score"))
    {
        opt.is_tier2_min_single_align_score=true;
    }

    if (vm.count("tier2-min-paired-align-score"))
    {
        opt.is_tier2_min_paired_align_score=true;
    }

    if (vm.count("tier2-mismatch-density-filter-count"))
    {
        opt.is_tier2_mismatch_density_filter_count=true;
    }

    check_option_arg_range(opt.tier2_indel_nonsite_match_prob,"tier2-indel-nonsite-match-prob",0.,1.,errorMsg);
    return (! errorMsg.empty());
}
