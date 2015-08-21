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


#include "starling_option_parser.hh"



po::options_description
get_starling_option_parser(
    starling_options& opt)
{
    po::options_description gvcf_opt("gVCF options");
    gvcf_opt.add_options()
    ("gvcf-file", po::value(&opt.gvcf.out_file),
     "gVCF output file-name, if not supplied output will be written to stdout.")
    ("chrom-depth-file", po::value(&opt.gvcf.chrom_depth_file),
     "If provided, the mean depth for each chromosome will be read from file, and these values will be used for high depth filtration. File should contain one line per chromosome, where each line begins with: \"chrom_name<TAB>depth\" (default: no chrom depth filtration)")
    ("gvcf-max-depth-factor", po::value(&opt.gvcf.max_depth_factor)->default_value(opt.gvcf.max_depth_factor),
     "If a chrom depth file is supplied then loci with depth exceeding the mean chromosome depth times this value are filtered")
    ("gvcf-min-gqx", po::value(&opt.gvcf.min_gqx)->default_value(opt.gvcf.min_gqx),
     "Minimum locus GQX in gVCF output. Providing a negative value disables the filter.")
    ("gvcf-max-snv-strand-bias", po::value(&opt.gvcf.max_snv_sb)->default_value(opt.gvcf.max_snv_sb),
     "Maximum SNV strand bias value")
    ("gvcf-no-snv-strand-bias-filter", po::value(&opt.gvcf.is_max_snv_sb)->zero_tokens()->implicit_value(false),
     "Disable SNV strand-bias filter")
    ("gvcf-max-snv-hpol", po::value(&opt.gvcf.max_snv_hpol)->default_value(opt.gvcf.max_snv_hpol),
     "SNVs are filtered if they exist in a homopolymer context greater than this length. A negative value disables the filter")
    ("gvcf-max-indel-ref-repeat", po::value(&opt.gvcf.max_ref_rep)->default_value(opt.gvcf.max_ref_rep),
     "Indels are filtered if they lengthen or contract a homopolymer or dinucleotide with reference repeat length greater than this value. A negative value disables the filter")
    ("gvcf-min-blockable-nonref", po::value(&opt.gvcf.block_max_nonref)->default_value(opt.gvcf.block_max_nonref),
     "A site cannot be joined into a non-variant block if it contains more than this fraction of non-reference alleles")
    ("gvcf-include-hapscore", po::value(&opt.is_compute_hapscore)->zero_tokens(),
     "Include haplotype score at SNV positions in gVCF output.")

    ("gvcf-block-percent-tol", po::value(&opt.gvcf.block_percent_tol)->default_value(opt.gvcf.block_percent_tol),
     "Non-variant blocks are chosen to constrain sample values to range [x,y], y <= max(x+3,x*(100+block-percent-tol)/100)")
    ("gvcf-no-block-compression", po::value(&opt.gvcf.is_block_compression)->zero_tokens()->implicit_value(false),
     "Turn off block compression in gVCF output")
    ("gvcf-report-VQSRmetrics", po::value(&opt.is_report_germline_VQSRmetrics)->zero_tokens(),
     "Report metrics used for germline VQSR")
    ("gvcf-compute-calibration-features", po::value(&opt.is_compute_calibration_features)->zero_tokens(),
     "Output all features used for calibration model training, development only.")
    ("nocompress-bed",  po::value(&opt.gvcf.nocompress_region_bedfile),
     "Bed file with sites that should not be block-compressed in gVCF (must be bgzip compressed and tabix indexed).")
    ("targeted-regions-bed",  po::value(&opt.gvcf.targeted_regions_bedfile),
     "Bed file with targeted regions. Variants outside these regions will be filtered. (must be bgzip compressed and tabix indexed).")
    ("indel-error-model",  po::value(&opt.indel_error_model_name)->default_value("new"),
     "Choose indel error model to use, available option old,new, new_stratified (development option only)")
    ("indel-ref-error-factor",  po::value(&opt.indel_ref_error_factor)->default_value(opt.indel_ref_error_factor),
     "Choose multiplier for ref error rate to use; 1 would be expected to be correct, but higher values counteract a bias away from homozygous indels (undercalling)")

    ("gvcf-skip-header", po::value(&opt.gvcf.is_skip_header)->zero_tokens(),
     "Skip writing header info for the gvcf file (usually used to simplify segment concatenation)")
    ;

    po::options_description phase_opt("Read-backed phasing options");
    phase_opt.add_options()
    ("do-short-range-phasing", po::value(&opt.do_codon_phasing)->zero_tokens(),
     "Enable short-range SNP phasing")
    ("phasing-window", po::value(&opt.phasing_window)->default_value(opt.phasing_window),
     "The maximum window to consider for short-range phasing")
    ;

    po::options_description score_opt("scoring-options");
    score_opt.add_options()
    ("variant-scoring-models-file", po::value(&opt.germline_variant_scoring_models_filename),
     "Model file for germline small variant scoring (VQSR)")
    ("variant-scoring-model-name", po::value(&opt.germline_variant_scoring_model_name),
     "The scoring model for germline small variants")
    ;

    po::options_description starling_parse_opt("Germline calling options");
    starling_parse_opt.add(gvcf_opt).add(phase_opt).add(score_opt);

    // final assembly
    po::options_description visible("Options");
    visible.add(starling_parse_opt);

    po::options_description visible2(get_starling_base_option_parser(opt));
    visible.add(visible2);

    po::options_description help_parse_opt("Help");
    help_parse_opt.add_options()
    ("help,h","print this message");

    visible.add(help_parse_opt);

    return visible;
}



void
finalize_starling_options(
    const prog_info& pinfo,
    const po::variables_map& vm,
    starling_options& opt)
{
    // gvcf option handlers:
    opt.gvcf.is_min_gqx = (opt.gvcf.min_gqx >= 0);
    opt.gvcf.is_max_snv_hpol = (opt.gvcf.max_snv_hpol >= 0);
    opt.gvcf.is_max_ref_rep = (opt.gvcf.max_ref_rep >= 0);

    if (opt.gvcf.block_percent_tol > 100)
    {
        pinfo.usage("block-percent-tol must be in range [0-100].");
    }

    finalize_starling_base_options(pinfo,vm,opt);
}

