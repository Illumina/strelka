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

/// \file
/// \author Chris Saunders
///


#include "starling_option_parser.hh"
#include "options/AlignmentFileOptionsParser.hh"

#include "boost/filesystem.hpp"


//#define DEBUG_OPTIONS
#ifdef DEBUG_OPTIONS
#include "blt_util/log.hh"
#endif




po::options_description
get_starling_option_parser(
    starling_options& opt)
{
    po::options_description aligndesc(getOptionsDescription(opt.alignFileOpt));

    po::options_description gvcf_opt("gVCF options");
    gvcf_opt.add_options()
    ("gvcf-output-prefix", po::value(&opt.gvcf.outputPrefix),
     "Turn on gvcf output mode and use the supplied prefix for all variant output files.")
    ("chrom-depth-file", po::value(&opt.gvcf.chrom_depth_file),
     "If provided, the mean depth for each chromosome will be read from file, and these values will be used for high depth filtration. File should contain one line per chromosome, where each line begins with: \"chrom_name<TAB>depth\" (default: no chrom depth filtration)")
    ("gvcf-max-depth-factor", po::value(&opt.gvcf.max_depth_factor)->default_value(opt.gvcf.max_depth_factor),
     "If a chrom depth file is supplied then loci with depth exceeding the mean chromosome depth times this value are filtered")
    ("gvcf-min-gqx", po::value(&opt.gvcf.min_gqx)->default_value(opt.gvcf.min_gqx),
     "Default minimum locus GQX in gVCF output (used for non-homref variants not passed through EVS). Providing a negative value disables the filter.")
    ("gvcf-min-homref-gqx", po::value(&opt.gvcf.min_homref_gqx)->default_value(opt.gvcf.min_homref_gqx),
     "Minimum homref locus GQX in gVCF output. Providing a negative value disables the filter.")
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

    ("gvcf-block-percent-tol", po::value(&opt.gvcf.block_percent_tol)->default_value(opt.gvcf.block_percent_tol),
     "Non-variant blocks are chosen to constrain sample values to range [x,y], y <= max(x+3,x*(100+block-percent-tol)/100)")
    ("gvcf-no-block-compression", po::value(&opt.gvcf.is_block_compression)->zero_tokens()->implicit_value(false),
     "Turn off block compression in gVCF output")
    ("nocompress-bed",  po::value(&opt.gvcf.nocompress_region_bedfile),
     "Bed file with sites that should not be block-compressed in gVCF (must be bgzip compressed and tabix indexed).")
    ("call-continuous-vf",  po::value(&opt.is_ploidy_prior)->zero_tokens()->implicit_value(false),
     "Instead of a haploid/diploid prior assumption, output a continuous VF")
    ("gvcf-skip-header", po::value(&opt.gvcf.is_skip_header)->zero_tokens(),
     "Skip writing header info for the gvcf file (usually used to simplify segment concatenation)")
    ("gvcf-include-header", po::value(&opt.gvcf.include_headers)->multitoken(),
     "Include the specified field description in the header (usually used to simplify segment concatenation when different segments have different fields)")
    ;

    po::options_description phase_opt("Read-backed phasing options");
    phase_opt.add_options()
    ("enable-read-backed-phasing", po::value(&opt.isUseVariantPhaser)->zero_tokens(),
     "Enable read-backed variant phasing")
    ;

    po::options_description score_opt("scoring-options");
    score_opt.add_options()
    ("snv-scoring-model-file", po::value(&opt.snv_scoring_model_filename),
     "File providing snv empirical scoring model")
    ("indel-scoring-model-file", po::value(&opt.indel_scoring_model_filename),
     "File providing indel empirical scoring model")
    ("use-rna-scoring", po::value(&opt.isRNA)->zero_tokens(),
     "Change to RNA-Seq analysis settings")
    ;

    po::options_description other_opt("other-options");
    other_opt.add_options()
    ("het-variant-frequency-extension", po::value(&opt.hetVariantFrequencyExtension)->default_value(opt.hetVariantFrequencyExtension),
     "Heterozygous variant allele frequency will be modeled as a range around 0.5 +/- this value. Value must be in [0,0.5)")
    ;

    po::options_description starling_parse_opt("Germline calling options");
    starling_parse_opt.add(aligndesc).add(gvcf_opt).add(phase_opt).add(score_opt).add(other_opt);

    // final assembly
    po::options_description visible("Options");
    visible.add(starling_parse_opt);

    po::options_description visible2(get_starling_base_option_parser(opt));
    visible.add(visible2);

    return visible;
}



void
finalize_starling_options(
    const prog_info& pinfo,
    const po::variables_map& vm,
    starling_options& opt)
{
    parseOptions(vm, opt.alignFileOpt);
    std::string errorMsg;
    if (checkOptions(opt.alignFileOpt, errorMsg))
    {
        pinfo.usage(errorMsg.c_str());
        //usage(log_os,prog,visible,errorMsg.c_str());
    }

    // gvcf option handlers:
    opt.gvcf.is_min_gqx = (opt.gvcf.min_gqx >= 0);
    opt.gvcf.is_min_homref_gqx = (opt.gvcf.min_homref_gqx >= 0);
    opt.gvcf.is_max_snv_hpol = (opt.gvcf.max_snv_hpol >= 0);

    if (opt.gvcf.block_percent_tol > 100)
    {
        pinfo.usage("block-percent-tol must be in range [0-100].");
    }

    if ((opt.hetVariantFrequencyExtension < 0.) || (opt.hetVariantFrequencyExtension >= 0.5))
    {
        pinfo.usage("het-variant-frequency-extension must be in range [0,0.5)\n");
    }

    if (opt.isReportEVSFeatures)
    {
        /// EVS feature output is contrained to the single-sample input case right now:
        const unsigned sampleCount(opt.alignFileOpt.alignmentFilenames.size());
        if (1 != sampleCount)
        {
            pinfo.usage("EVS features can only be reported when analyzing a single sample");
        }
    }

    checkOptionalInputFile(pinfo, opt.snv_scoring_model_filename, "SNV empirical scoring model");
    checkOptionalInputFile(pinfo, opt.indel_scoring_model_filename, "Indel empirical scoring model");

    finalize_starling_base_options(pinfo,vm,opt);
}
