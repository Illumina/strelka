// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
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


#include "starling_option_parser.hh"

//#define DEBUG_OPTIONS
#ifdef DEBUG_OPTIONS
#include "blt_util/log.hh"
#endif




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
    ("nocompress-bed",  po::value(&opt.gvcf.nocompress_region_bedfile),
     "Bed file with sites that should not be block-compressed in gVCF (must be bgzip compressed and tabix indexed).")
    ("targeted-regions-bed",  po::value(&opt.gvcf.targeted_regions_bedfile),
     "Bed file with targeted regions. Variants outside these regions will be filtered. (must be bgzip compressed and tabix indexed).")
    ("indel-error-model",  po::value(&opt.indel_error_model_name)->default_value("new"),
     "Choose indel error model to use, available option old,new, new_stratified (development option only)")
    ("indel-ref-error-factor",  po::value(&opt.indel_ref_error_factor)->default_value(opt.indel_ref_error_factor),
     "Choose multiplier for ref error rate to use; 1 would be expected to be correct, but higher values counteract a bias away from homozygous indels (undercalling)")
    ("call-continuous-vf",  po::value(&opt.is_ploidy_prior)->zero_tokens()->implicit_value(false),
     "Instead of a haploid/diploid prior assumption, output a continuous VF")
    ("noise-floor",  po::value(&opt.noise_floor)->default_value(opt.noise_floor),
     "The noise rate for basecalls assumed when calling continuous variant frequencies")
    ("min-het-vf",  po::value(&opt.min_het_vf)->default_value(opt.min_het_vf),
     "The minimum allele frequency to call a heterozygous genotype when calling continuous variant frequencies")
    ("gvcf-skip-header", po::value(&opt.gvcf.is_skip_header)->zero_tokens(),
     "Skip writing header info for the gvcf file (usually used to simplify segment concatenation)")
    ("gvcf-include-header", po::value(&opt.gvcf.include_headers)->multitoken(),
     "Include the specified field description in the header (usually used to simplify segment concatenation when different segments have different fields)")
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
     "File providing germline empirical variant scoring models")
    ("variant-scoring-model-name", po::value(&opt.germline_variant_scoring_model_name),
     "The scoring model for germline variants")
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
    if (!opt.is_ploidy_prior)
    {
        if (opt.noise_floor < 0)
        {
            // not specified - derive from the min qscore
            opt.noise_floor = pow(10, (-1.0 * opt.min_qscore)/10.0);
#ifdef DEBUG_OPTIONS
            log_os << "Setting noise floor to: " << opt.noise_floor << " from min_qscore " << opt.min_qscore << "\n";
#endif
        }
        else if (opt.noise_floor > 0 && opt.noise_floor < 0.5)
        {
            // specified explicitly. Override the min_qscore if it is incompatible with the noise floor
            int min_qscore = (int)round(-10 * log10(opt.noise_floor));
            if (min_qscore > opt.min_qscore)
            {
                opt.min_qscore = min_qscore;
#ifdef DEBUG_OPTIONS
                log_os << "Overriding min q_score to: " << min_qscore << " to support noise floor " << opt.noise_floor << "\n";
#endif
            }

        }
        if (opt.noise_floor >= 0.5 || opt.noise_floor <= 0.0)
        {
            pinfo.usage("noise-floor must be in range (0, 0.5)");
        }
        if (opt.min_het_vf <= 0.0 || opt.min_het_vf >= 0.5)
        {
            pinfo.usage("min-het-vf must be in range (0, 0.5)");
        }

    }

    finalize_starling_base_options(pinfo,vm,opt);
}

