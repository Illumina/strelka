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

#include "strelka_option_parser.hh"

#include "blt_common/blt_arg_validate.hh"
#include "options/AlignmentFileOptionsParser.hh"
#include "options/TumorNormalAlignmentFileOptionsParser.hh"
#include "starling_common/starling_base_option_parser.hh"
#include "starling_common/Tier2OptionsParser.hh"



po::options_description
get_strelka_option_parser(
    strelka_options& opt)
{
    po::options_description aligndesc(getOptionsDescription(opt.alignFileOpt));

    po::options_description strelka_parse_opt_sv("Somatic variant-calling");
    strelka_parse_opt_sv.add_options()
    ("somatic-snv-file",
     po::value(&opt.somatic_snv_filename),
     "Output file for somatic snv-calls (note this uses settings from the bsnp diploid caller for the normal sample)")
    ("somatic-snv-rate",
     po::value(&opt.somatic_snv_rate)->default_value(opt.somatic_snv_rate),
     "Expected rate of somatic snvs (allowed range: [0-1])")
    ("somatic-indel-file",
     po::value(&opt.somatic_indel_filename),
     "Output file for somatic indel (note this uses settings from the bindel diploid caller for the normal sample)")
    ("somatic-indel-rate",
     po::value(&opt.somatic_indel_rate)->default_value(opt.somatic_indel_rate),
     "Expected rate of somatic indels (allowed range: [0-1])")
    ("shared-site-error-rate",
     po::value(&opt.shared_site_error_rate)->default_value(opt.shared_site_error_rate),
     "Expected rate of site specific errors shared in the tumor and normal data.")
    ("shared-indel-error-factor",
     po::value(&opt.shared_indel_error_factor)->default_value(opt.shared_indel_error_factor),
     "Factor affecting the expected rate of context-specific spurious indel errors shared in the tumor and normal data.")
    ("shared-site-error-strand-bias-fraction",
     po::value(&opt.shared_site_error_strand_bias_fraction)->default_value(opt.shared_site_error_strand_bias_fraction),
     "Expected fraction of site-specific errors which are single-stranded.")
    ("ssnv-contam-tolerance",
     po::value(&opt.ssnv_contam_tolerance)->default_value(opt.ssnv_contam_tolerance),
     "Tolerance of tumor contamination in the normal sample for SNVs (allowed range: [0-1]).")
    ("indel-contam-tolerance",
     po::value(&opt.indel_contam_tolerance)->default_value(opt.indel_contam_tolerance),
     "Tolerance of tumor contamination in the normal sample for indels (allowed range: [0-1]).")
    ("somatic-callable-regions-file", po::value(&opt.somatic_callable_filename),
     "Output a bed file of regions which are confidently somatic or non-somatic for SNVs at allele frequencies of 10% or greater.")
    ("noise-vcf", po::value(&opt.noise_vcf)->multitoken(),
     "Noise panel VCF for low-frequency noise")
    ;

    po::options_description strelka_parse_opt_filter("Somatic variant-calling filters");
    strelka_parse_opt_filter.add_options()
    ("strelka-chrom-depth-file", po::value(&opt.sfilter.chrom_depth_file),
     "If provided, the expected depth for each chromosome will be read from file, and these values will be used for high depth filtration. File should contain one line per chromosome, where each line begins with: \"chrom_name<TAB>depth\" (default: no chrom depth filtration)")
    ("strelka-max-depth-factor", po::value(&opt.sfilter.max_depth_factor)->default_value(opt.sfilter.max_depth_factor),
     "If a chrom depth file is supplied then loci with depth exceeding the mean chromosome depth times this value are filtered")
    ("strelka-skip-header", po::value(&opt.sfilter.is_skip_header)->zero_tokens(),
     "Skip writing header info for the somatic vcf/bed files (usually used to simplify segment concatenation)")
    // snv only:
    ("strelka-snv-max-filtered-basecall-frac", po::value(&opt.sfilter.snv_max_filtered_basecall_frac)->default_value(opt.sfilter.snv_max_filtered_basecall_frac),
     "max filtered call fraction")
    ("strelka-snv-max-spanning-deletion-frac", po::value(&opt.sfilter.snv_max_spanning_deletion_frac)->default_value(opt.sfilter.snv_max_spanning_deletion_frac),
     "max fraction of overlapping deletion reads")
    ("strelka-snv-min-qss-ref", po::value(&opt.sfilter.snv_min_qss_ref)->default_value(opt.sfilter.snv_min_qss_ref),
     "min QSS_ref value")
    // indel only:
    ("strelka-indel-max-window-filtered-basecall-frac",  po::value(&opt.sfilter.indelMaxWindowFilteredBasecallFrac)->default_value(opt.sfilter.indelMaxWindowFilteredBasecallFrac),
     "indel are filtered if more than this fraction of basecalls are filtered in a 50 base window")
    ("strelka-indel-min-qsi-ref", po::value(&opt.sfilter.sindelQuality_LowerBound)->default_value(opt.sfilter.sindelQuality_LowerBound),
     "min QSI_ref value")
    ;

    po::options_description tier2_opt(getTier2OptionsDescription(opt.tier2));

    po::options_description score_opt("scoring-options");
    score_opt.add_options()
    ("somatic-snv-scoring-model-file", po::value(&opt.somatic_snv_scoring_model_filename),
     "Model file for somatic SNV scoring")
    ("somatic-indel-scoring-model-file", po::value(&opt.somatic_indel_scoring_model_filename),
     "Model file for somatic indel scoring")
    ;

    po::options_description strelka_parse_opt("Two-sample options");
    strelka_parse_opt
    .add(aligndesc).add(strelka_parse_opt_sv).add(strelka_parse_opt_filter)
    .add(tier2_opt).add(score_opt);

    // final assembly
    po::options_description visible("Options");
    visible.add(strelka_parse_opt);

    // add starling base options:
    po::options_description visible2(get_starling_base_option_parser(opt));
    visible.add(visible2);

    return visible;
}



void
finalize_strelka_options(
    const prog_info& pinfo,
    const po::variables_map& vm,
    strelka_options& opt)
{
    parseOptions(vm, opt.alignFileOpt);
    std::string errorMsg;
    if (checkOptions(opt.alignFileOpt, errorMsg))
    {
        pinfo.usage(errorMsg.c_str());
        //usage(log_os,prog,visible,errorMsg.c_str());
    }

    {
        unsigned normalCount(0);
        unsigned tumorCount(0);
        for (const bool isTumor : opt.alignFileOpt.isAlignmentTumor)
        {
            if (isTumor)
            {
                tumorCount++;
            }
            else
            {
                normalCount++;
            }
        }

        // undocumented research option to provide zero normal input. somatic caller will run but no sensible scores yet
        if (normalCount > 1)
        {
            pinfo.usage("Must specify no more than one normal sample alignment file.");
        }
        if (tumorCount != 1)
        {
            pinfo.usage("Must specify exactly one tumor sample alignment file.");
        }
    }

    check_option_arg_range(pinfo,opt.somatic_snv_rate,"somatic-snv-rate",0.,1.);
    check_option_arg_range(pinfo,opt.shared_site_error_rate,"shared-site-error-rate",0.,1.);
    check_option_arg_range(pinfo,opt.shared_site_error_strand_bias_fraction,"shared-site-strand-strand-bias-fraction",0.,1.);

    check_option_arg_range(pinfo,opt.somatic_indel_rate,"somatic-indel-rate",0.,1.);

    // deal with sfilter options:
    if (opt.sfilter.max_depth_factor < 0)
    {
        pinfo.usage("Strelka depth factor must not be less than 0");
    }

    checkOptionalInputFile(pinfo, opt.somatic_snv_scoring_model_filename, "somatic snv scoring model");
    checkOptionalInputFile(pinfo, opt.somatic_indel_scoring_model_filename, "somatic indel scoring model");

    finalize_starling_base_options(pinfo,vm,opt);
}
