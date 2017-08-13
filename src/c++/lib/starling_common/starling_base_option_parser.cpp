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

#include "starling_base_option_parser.hh"

#include "blt_common/blt_arg_validate.hh"
#include "blt_util/compat_util.hh"
#include "starling_common/Tier2OptionsParser.hh"

#include "boost/filesystem.hpp"
#include "boost/format.hpp"

#include <iostream>
#include <sstream>

//#define DEBUG_OPTIONPARSER

#ifdef DEBUG_OPTIONPARSER
#include "blt_util/log.hh"
#endif



void
checkOptionalFile(
    const prog_info& pinfo,
    const std::string& filename,
    const char* label)
{
    if (filename.empty()) return;
    if (! boost::filesystem::exists(filename))
    {
        std::ostringstream oss;
        oss << "Can't find " << label << " file '" << filename << "'";
        pinfo.usage(oss.str().c_str());
    }
}



po::options_description
get_starling_base_option_parser(
    starling_base_options& opt)
{
    po::options_description core_opt("core options");
    core_opt.add_options()
    ("ref", po::value(&opt.referenceFilename),
     "fasta reference sequence, samtools index file must be present (required)")
    ("region", po::value<regions_t>(),
     "samtools formatted region, eg. 'chr1:20-30'. May be supplied more than once but regions must not overlap. At least one entry required.")
    ;

    po::options_description geno_opt("genotyping options");
    geno_opt.add_options()
    ("snp-theta", po::value(&opt.bsnp_diploid_theta)->default_value(opt.bsnp_diploid_theta),
     "Set snp theta.")
    ("indel-theta", po::value(&opt.bindel_diploid_theta)->default_value(opt.bindel_diploid_theta),
     "Set indel theta")
    ;

    po::options_description realign_opt("realignment-options");
    realign_opt.add_options()
    ("max-indel-toggle-depth", po::value(&opt.max_read_indel_toggle)->default_value(opt.max_read_indel_toggle),
     "Controls the realignment stringency. Lowering this value will increase the realignment speed at the expense of indel-call quality")
    ("retain-optimal-soft-clipping", po::value(&opt.isRetainOptimalSoftClipping)->zero_tokens(),
     "Retain input alignment soft-clipping if it outscores realignment with soft-clipping unrolled.")
    ;

    po::options_description indel_opt("indel-options");
    indel_opt.add_options()
    ("max-candidate-indel-depth",
     po::value(&opt.max_candidate_indel_depth)->default_value(opt.max_candidate_indel_depth),
     "Maximum estimated read depth for an indel to reach candidacy. If any one sample exceeds this depth at the indel, the indel will not reach candidacy in all indel-syncronized samples. A non-positive value disables the filter")
    ("max-candidate-indel-depth-factor",
     po::value(&opt.max_candidate_indel_depth_factor)->default_value(opt.max_candidate_indel_depth_factor),
     "If a chromosome maximum depth filter is in use, then at this factor of the filtration depth cutoff no indels will reach candidacy in all indel-synchronized samples. A non-positive value disables the filter")
    ("min-candidate-open-length",
     po::value(&opt.min_candidate_indel_open_length)->default_value(opt.min_candidate_indel_open_length),
     "Minimum open-ended breakpoint sequence length required to become a breakpoint candidate")
    ("candidate-indel-input-vcf",
     po::value(&opt.input_candidate_indel_vcf)->multitoken(),
     "Add candidate indels from the specified vcf file. Option can be provided multiple times to combine evidence from multiple vcf files. Variants will be ignored if not correctly normalized. (must be bgzip compressed and tabix indexed)")
    ("force-output-vcf", po::value(&opt.force_output_vcf)->multitoken(),
     "Force each site or indel in the vcf file to be written to the snv or indel output, even if no variant is found. Any indels submitted will also be treated as candidate indels. Option can be provided multiple times to combine multiple vcf files. Unnormalized variants will trigger a runtime error. (must be bgzip compressed and tabix indexed)")
    ("upstream-oligo-size", po::value(&opt.upstream_oligo_size),
     "Treat reads as if they have an upstream oligo anchor for purposes of meeting minimum breakpoint overlap in support of an indel.")
    ;

    po::options_description ploidy_opt("ploidy-options");
    ploidy_opt.add_options()
    ("ploidy-region-vcf", po::value(&opt.ploidy_region_vcf),
     "Specify vcf file describing ploidy of regions. Regions span [POS+1,INFO/END] and use FORMAT/CN to indicate the per-sample region ploidy. Any CN value besides 1 and 0 are ignored. (must be bgzip compressed and tabix indexed)")
    ;

    po::options_description input_opt("input-options");
    input_opt.add_options()
    ("max-input-depth", po::value(&opt.max_input_depth),
     "Maximum allowed read depth per sample (prior to realignment). Input reads which would exceed this depth are filtered out.  (default: no limit)")
    ("max-sample-read-buffer", po::value(&opt.maxBufferedReads)->default_value(opt.maxBufferedReads),
     "Maximum reads buffered for each sample")
    ;

    po::options_description other_opt("other-options");
    other_opt.add_options()
    ("stats-file", po::value(&opt.segmentStatsFilename),
     "Write runtime stats to file")
    ("report-evs-features", po::value(&opt.isReportEVSFeatures)->zero_tokens(),
     "Report empirical variant scoring (EVS) training features in VCF output")
    ("indel-error-models-file", po::value<std::vector<std::string>>(&opt.indelErrorModelFilenames),
     "File containing indel error models (may be specified more than once)")
    ("theta-file", po::value<std::string>(&opt.thetaFilename),
     "File containing theta values")
    ("indel-error-model-name", po::value(&opt.indel_error_model_name)->default_value(opt.indel_error_model_name),
     "Static indel error model name. If no indel error model file is provided a hard-coded model can be selected with this argument instead. Current options are ('adaptiveDefault','logLinear'). This option is ignored when at least one indel error models file is provided.")
    ("call-regions-bed",  po::value(&opt.callRegionsBedFilename),
     "Bed file describing regions to call. No output will be provided outside of these regions. (must be bgzip compressed and tabix indexed).")
    ;

    po::options_description new_opt("Shared small-variant options");

    new_opt.add(core_opt).add(geno_opt);
    new_opt.add(realign_opt).add(indel_opt).add(ploidy_opt);
    new_opt.add(input_opt).add(other_opt);

    return new_opt;
}



void
write_starling_legacy_options(
    const starling_base_options& default_opt,
    std::ostream& os)
{
    os <<
       " -bsnp-diploid-het-bias x\n"
       "                    - Set bias term for the heterozygous state in the bsnp model, such that\n"
       "                      hets are expected at allele ratios in the range [0.5-x,0.5+x] (default: 0)\n"
#if 0
       " -bsnp-nploid n x   - Use Bayesian nploid genotype snp caller with ploidy=n and prior(snp)=x\n"
       " -lsnp-alpha x      - Use likelihood ratio test snp caller with alpha=x\n"
#endif
       " -min-qscore n      - Don't use base if qscore<n (default: " << default_opt.min_qscore << ")\n"
       " -max-window-mismatch n m\n"
       "                    - Don't use base if mismatch count>n within a window of m flanking bases\n"
       " -min-mapping-quality n\n"
       "                    - Reads with mapping quality<n are marked as tier1 filtered. Such reads are not\n"
       "                      directly used for variant calling unless a tier2 value is defined in certain applications.\n"
       "                      Filtered reads may also still be used to compute locus quality metrics (default score: " << default_opt.min_mapping_quality << ")\n"
       " -include-singleton - Include paired reads with unmapped mates\n"
       " -include-anomalous - Include paired reads which are not part of a 'proper pair' (anomalous orientation or insert size)\n"
       " -print-evidence    - Print the observed data at single site events (does not include indels)\n"
       " -print-all-site-evidence\n"
       "                    - Print the observed data for all sites (does not include indels)\n"
       " -bindel-diploid-file file\n"
       "                    - Run Bayesian diploid genotype caller, write results to 'file'\n"
       " -indel-error-rate x\n"
       "                    - If calling indels, set the indel error rate to a constant value of x (0<=x<=1).\n"
       "                      The default indel error rate is taken from an empirical function accounting for\n"
       "                      homopolymer length and indel type (i.e. insertion or deletion). This option\n"
       "                      overrides the default behavior.\n"
       " -indel-nonsite-match-prob x\n"
       "                    - The probability of a base matching the reference in an 'average' mismapped read. This\n"
       "                      value is used by the indel-caller only. (default: " << default_opt.indel_nonsite_match_prob << ")\n"
       " -genome-size n     - Specify the total number of non-ambiguous bases in the genome to which the input reads\n"
       "                      have been aligned for use in indel calling.\n"
       " -max-candidate-indel-density x\n"
       "                    - If there are more than x candidate indels per base intersecting a read, then realignment\n"
       "                      is truncated to only allow individual indel toggles of the starting alignments for that read.\n"
       "                      x must be greater than 0  (default: " << default_opt.max_candidate_indel_density << ")\n"
       " -realigned-read-file file\n"
       "                    - Write reads which have had their alignments altered during realignemnt to a BAM file.\n"
       " -realign-submapped-reads\n"
       "                    - When this argument is provided, even reads which fail the variant calling mapping thresholds\n"
       "                      are realigned using the same procedure as the variant calling reads.\n"
       " -no-ambiguous-path-clip\n"
       "                    - Turn off ambiguous read trimming after realignment.\n"
       " -max-indel-size    - Sets the maximum size for indels processed for indel genotype calling and realignment.\n"
       "                      Increasing this value should lead to an approx linear increase in memory consumption.\n"
       "                      (default: " << default_opt.maxIndelSize << ")\n"
       "\n"
       " -all-warnings      - print all warnings (default: errors and low-frequency warnings only)\n"
       "\n"
       " -h                 - Display usage (this page)\n";
}



static
void
finalize_legacy_starling_options(
    const prog_info& pinfo,
    starling_base_options& opt)
{
    if (! opt.is_user_genome_size)
    {
        // this requirement is not what we want, but it's the only way to make things reliable for now:
        pinfo.usage("must specify genome-size");
    }
    else
    {
        if (opt.user_genome_size<1)
        {
            pinfo.usage("genome-size must be greater than 0");
        }
    }

    validate_blt_opt(pinfo,opt);
}



void
finalize_starling_base_options(
    const prog_info& pinfo,
    const po::variables_map& vm,
    starling_base_options& opt)
{
    // check for and sanitize the reference fasta sequence
    if (opt.referenceFilename.empty())
    {
        pinfo.usage("Must specify a fasta reference file");
    }

    // canonicalize the reference sequence path:
    /// TODO: replace this with the same thing from boost?
    if (! compat_realpath(opt.referenceFilename))
    {
        std::ostringstream oss;
        oss << "can't resolve reference path: " << opt.referenceFilename << "\n";
        pinfo.usage(oss.str().c_str());
    }

    // set analysis regions:
    if (vm.count("region"))
    {
        opt.regions=(boost::any_cast<regions_t>(vm["region"].value()));
    }

    // validate regions:
    {
        if (opt.regions.empty())
        {
            pinfo.usage("Need at least one samtools formatted region");
        }
        for (const auto& region : opt.regions)
        {
            if (region.empty())
            {
                pinfo.usage("Empty region argument");
            }
        }
    }

    // max_theta for indels is actually 2./3., but because we don't
    // allow non-reference hets, we stick with the lower value
    // used for snps:
    //
    if (opt.bindel_diploid_theta>MAX_DIPLOID_THETA)
    {
        std::ostringstream oss;
        oss << "indel diploid heterozygosity exceeds maximum value of: " << MAX_DIPLOID_THETA;
        pinfo.usage(oss.str().c_str());
    }

    if (vm.count("max-input-depth"))
    {
        opt.is_max_input_depth=true;
    }

    for (const auto& indelErrorModelFilename : opt.indelErrorModelFilenames)
    {
        checkOptionalFile(pinfo, indelErrorModelFilename, "indel error models");
    }
    /// tier2 options are not parsed by starling_base, but need to live up here for now,
    /// so validate them together with the rest of starling_base
    std::string errorMsg;
    if (parseTier2Options(vm,opt.tier2,errorMsg))
    {
        pinfo.usage(errorMsg.c_str());
    }

    if (opt.tier2.is_tier2_min_mapping_quality)
    {
        if (opt.tier2.tier2_min_mapping_quality >= opt.min_mapping_quality)
        {
            std::ostringstream oss;
            oss << "Invalid tier2 min mapping quality. Value must be lower than tier1 min mapping quality: '" << opt.min_mapping_quality << "'";
            pinfo.usage(oss.str().c_str());
        }
    }

    finalize_legacy_starling_options(pinfo,opt);
}
