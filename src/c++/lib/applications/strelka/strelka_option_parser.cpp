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

#include "strelka_option_parser.hh"

#include "blt_common/blt_arg_validate.hh"
#include "starling_common/starling_option_parser.hh"



po::options_description
get_strelka_option_parser(
    strelka_options& opt)
{
    po::options_description strelka_parse_opt_ti("Tumor-sample input");
    strelka_parse_opt_ti.add_options()
    ("tumor-bam-file",
     po::value(&opt.tumor_bam_filename),
     "BAM/CRAM file containing read alignments for the tumor sample (required)")
    ;

    po::options_description strelka_parse_opt_to("Tumor-sample output");
    strelka_parse_opt_to.add_options()
    ("tumor-realigned-read-file",
     po::value(&opt.tumor_realigned_read_filename),
     "Write tumor reads which have had their alignments altered during realignment to a BAM file.")
    ;

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
#if 0
    ("somatic-indel-depth-window-file",
     po::value(&opt.somatic_indel_depth_window_filename),
     "Output file for depth window averages corresponding to each somatic indel call.")
    ("somatic-indel-depth-window-flank",
     po::value(&opt.depth_window_flank),
     "Somatic indel depth window flank-size (must be >1)")
#endif
    ("somatic-indel-rate",
     po::value(&opt.somatic_indel_rate)->default_value(opt.somatic_indel_rate),
     "Expected rate of somatic indels (allowed range: [0-1])")
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
    ("shared-site-error-rate",
     po::value(&opt.shared_site_error_rate)->default_value(opt.shared_site_error_rate),
     "Expected rate of site specific errors shared in the tumor and normal data.")
    ("shared-site-error-strand-bias-fraction",
     po::value(&opt.shared_site_error_strand_bias_fraction)->default_value(opt.shared_site_error_strand_bias_fraction),
     "Expected fraction of site-specific errors which are single-stranded.")
    ("site-somatic-normal-noise-rate",
     po::value(&opt.site_somatic_normal_noise_rate),
     "Expected rate of 'noise' in the normal sample at somatic call sites -- this allows for some degree of tumor contamination in the normal for raw somatic Q-scores (default: use shared site error instead)")
    ("shared-indel-error-rate",
     po::value(&opt.shared_indel_error_rate)->default_value(opt.shared_indel_error_rate),
     "Expected rate of site-specific spurious indel errors shared in the tumor and normal data.")
    ("tumor-min-candidate-indel-reads",
     po::value(&opt.tumor_sample_min_candidate_indel_reads),
     "Unless an indel is supported by at least this many reads in the tumor sample, it cannot become a candidate unless the global read count test passes for all samples. (default: not used)")
    ("tumor-min-small-candidate-indel-read-frac",
     po::value(&opt.tumor_sample_min_small_candidate_indel_read_frac),
     "For small indels an additional indel candidacy filter is applied: Unless at least this fraction of intersecting reads contain the small indel in the tumor sample, it cannot become a candidate unless this same test passes for other samples. (default: not used)")
    ("indel-somatic-normal-noise-rate",
     po::value(&opt.indel_somatic_normal_noise_rate),
     "Expected rate of 'noise' in the normal sample at somatic indels -- this allows for some degree of tumor contamination in the normal sample for raw somatic Q-scores (default: use shared site error instead)")
    ("somatic-callable-region-file",
     po::value(&opt.somatic_callable_filename),
     "Output a bed file of regions which are confidently somatic or non-somatic for SNVs at allele frequencies of 10% or greater.")
    ("noise-vcf", po::value(&opt.noise_vcf)->multitoken(),
     "Noise panel VCF for low-frequency noise")
    ;

    po::options_description strelka_parse_opt_filter("Somatic variant-calling filters");
    strelka_parse_opt_filter.add_options()
    ("strelka-chrom-depth-file", po::value(&opt.sfilter.chrom_depth_file),
     "If provided, the mean depth for each chromosome will be read from file, and these values will be used for high depth filtration. File should contain one line per chromosome, where each line begins with: \"chrom_name<TAB>depth\" (default: no chrom depth filtration)")
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
//     ("compute-VQSR-metrics", po::value(&opt.sfilter.compute_VQSR_options)->default_value(opt.sfilter.compute_VQSR_options),
//      "report VQSR metrics in variant records")
    // indel only:
    ("strelka-indel-max-ref-repeat",  po::value(&opt.sfilter.indelMaxRefRepeat)->default_value(opt.sfilter.indelMaxRefRepeat),
     "indels expand/contracting pattterns with greater than this repeat will be filtered out")
    ("strelka-indel-max-int-hpol-length",  po::value(&opt.sfilter.indelMaxIntHpolLength)->default_value(opt.sfilter.indelMaxIntHpolLength),
     "indels are filtered if they overlap 'interrupted homopolymers' greater than this length")
    ("strelka-indel-max-window-filtered-basecall-frac",  po::value(&opt.sfilter.indelMaxWindowFilteredBasecallFrac)->default_value(opt.sfilter.indelMaxWindowFilteredBasecallFrac),
     "indel are filtered if more than this fraction of basecalls are filtered in a 50 base window")
    ("strelka-indel-min-qsi-ref", po::value(&opt.sfilter.sindelQuality_LowerBound)->default_value(opt.sfilter.sindelQuality_LowerBound),
     "min QSI_ref value")
    ;

    po::options_description strelka_parse_opt("Two-sample options");
    strelka_parse_opt.add(strelka_parse_opt_ti).add(strelka_parse_opt_to).add(strelka_parse_opt_sv).add(strelka_parse_opt_filter);

    po::options_description help_parse_opt("Help");
    help_parse_opt.add_options()
    ("help,h","print this message");

    po::options_description visible("Options");
    visible.add(strelka_parse_opt).add(help_parse_opt);

    return visible;
}



void
finalize_strelka_options(const prog_info& pinfo,
                         const po::variables_map& vm,
                         strelka_options& opt)
{
    if (opt.tumor_bam_filename.empty())
    {
        pinfo.usage("Must specify a sorted & indexed BAM/CRAM file containing aligned tumor sample reads");
    }

    if (vm.count("skip-realignment"))
    {
        if (opt.is_call_indels())
        {
            pinfo.usage("Cannot disable realignment when indel-calling is selected.");
        }
    }

    check_option_arg_range(pinfo,opt.somatic_snv_rate,"somatic-snv-rate",0.,1.);
    check_option_arg_range(pinfo,opt.shared_site_error_rate,"shared-site-error-rate",0.,1.);
    check_option_arg_range(pinfo,opt.shared_site_error_strand_bias_fraction,"shared-site-strand-strand-bias-fraction",0.,1.);
    check_option_arg_range(pinfo,opt.shared_site_error_strand_bias_fraction,"site-somatic-normal-noise-rate",0.,1.);

    check_option_arg_range(pinfo,opt.somatic_indel_rate,"somatic-indel-rate",0.,1.);
    check_option_arg_range(pinfo,opt.shared_indel_error_rate,"shared-indel-error-rate",0.,1.);
    check_option_arg_range(pinfo,opt.shared_site_error_strand_bias_fraction,"indel-somatic-normal-noise-rate",0.,1.);
    check_option_arg_range(pinfo,opt.tier2_indel_nonsite_match_prob,"tier2-indel-nonsite-match-prob",0.,1.);

    if (vm.count("site-somatic-normal-noise-rate"))
    {
        opt.is_site_somatic_normal_noise_rate=true;
    }

    if (vm.count("indel-somatic-normal-noise-rate"))
    {
        opt.is_indel_somatic_normal_noise_rate=true;
    }

    if (vm.count("tier2-indel-nonsite-match-prob"))
    {
        opt.is_tier2_indel_nonsite_match_prob=true;
    }

    if (vm.count("tier2-min-single-align-score"))
    {
        if (opt.tier2_min_single_align_score >= opt.min_single_align_score)
        {
            std::ostringstream oss;
            oss << "Invalid tier2 single align score. Value must be lower than primary single align score: '" << opt.min_single_align_score << "'";
            pinfo.usage(oss.str().c_str());
        }
        opt.is_tier2_min_single_align_score=true;
    }

    if (vm.count("tier2-min-paired-align-score"))
    {
        if (opt.tier2_min_paired_align_score >= opt.min_paired_align_score)
        {
            std::ostringstream oss;
            oss << "Invalid tier2 paired align score. Value must be lower than primary paired align score: '" << opt.min_paired_align_score << "'";
            pinfo.usage(oss.str().c_str());
        }
        opt.is_tier2_min_paired_align_score=true;
    }

    if (vm.count("tier2-mismatch-density-filter-count"))
    {
        if (opt.is_tier2_no_mismatch_density_filter)
        {
            pinfo.usage("Only one tier2 mismatch density filter setting allowed");
        }
        opt.is_tier2_mismatch_density_filter_count=true;
    }

    if (vm.count("tumor-min-candidate-indel-reads"))
    {
        opt.is_tumor_sample_min_candidate_indel_reads=true;
    }

    if (vm.count("tumor-min-small-candidate-indel-read-frac"))
    {
        opt.is_tumor_sample_min_small_candidate_indel_read_frac=true;
    }

    // deal with sfilter options:
    if (opt.sfilter.max_depth_factor < 0)
    {
        pinfo.usage("Strelka depth factor must not be less than 0");
    }

    finalize_starling_options(pinfo,vm,opt);
}
