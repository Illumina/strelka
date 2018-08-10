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

#include "starling_base_option_parser.hh"

#include "blt_common/blt_arg_validate.hh"
#include "blt_util/compat_util.hh"
#include "options/optionsUtil.hh"
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
checkOptionalInputFile(
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
    ("realigned-output-prefix", po::value(&opt.realignedReadFilenamePrefix),
     "Write reads which have had their alignments altered during realignemnt to a BAM file, with the given path prefix.")
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
    ("max-indel-size", po::value(&opt.maxIndelSize)->default_value(opt.maxIndelSize),
     "Maximum size of indels processed for realignment and calling.")
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
    ("min-qscore", po::value(&opt.minBasecallErrorPhredProb)->default_value(opt.minBasecallErrorPhredProb),
     "Don't use a basecall for SNV calling if qscore is below this value.")
    ("min-mapping-quality", po::value(&opt.minMappingErrorPhredProb)->default_value(opt.minMappingErrorPhredProb),
     "Reads with mapping quality<n are marked as tier1 filtered. Such reads are not directly used for variant calling unless a tier2 value is defined in certain applications. Filtered reads may also still be used to compute locus quality metrics.")
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

    po::options_description help_parse_opt("Help");
    help_parse_opt.add_options()
    ("help,h","print this message");

    new_opt.add(help_parse_opt);

    return new_opt;
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

    // Convert reference sequence path to an absolute path. Just like for BAM/CRAM files, we want absolute but not
    // canonical path in this case, because following softlinks to the canonical path will often cause us to miss the
    // sidecar index file.
    {
        std::string errorMsg;
        if (checkAndStandardizeRequiredInputFilePath(opt.referenceFilename, "reference fasta", errorMsg))
        {
            pinfo.usage(errorMsg.c_str());
        }
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

    if (opt.bsnp_diploid_theta>MAX_DIPLOID_THETA)
    {
        std::ostringstream oss;
        oss << "diploid heterozygosity exceeds maximum value of: " << MAX_DIPLOID_THETA;
        pinfo.usage(oss.str().c_str());
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
        checkOptionalInputFile(pinfo, indelErrorModelFilename, "indel error models");
    }

    // tier2 options are not parsed by starling_base, but need to live up here for now,
    // so validate them together with the rest of starling_base
    {
        std::string errorMsg;
        if (parseTier2Options(vm, opt.tier2, errorMsg))
        {
            pinfo.usage(errorMsg.c_str());
        }
    }

    if (opt.useTier2Evidence)
    {
        if (opt.tier2.minMappingErrorPhredProb >= opt.minMappingErrorPhredProb)
        {
            std::ostringstream oss;
            oss << "Invalid tier2 min mapping quality. Value must be lower than tier1 min mapping quality: '" << opt.minMappingErrorPhredProb << "'";
            pinfo.usage(oss.str().c_str());
        }
    }
}
