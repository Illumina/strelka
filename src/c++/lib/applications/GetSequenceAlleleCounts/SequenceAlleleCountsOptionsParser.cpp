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

#include "SequenceAlleleCountsOptionsParser.hh"
#include "options/AlignmentFileOptionsParser.hh"


#include "boost/filesystem.hpp"



//#define DEBUG_OPTIONS
#ifdef DEBUG_OPTIONS
#include "blt_util/log.hh"
#endif




po::options_description
getSequenceAlleleCountsOptionsParser(
    SequenceAlleleCountsOptions& opt)
{
    po::options_description aligndesc(getOptionsDescription(opt.alignFileOpt));

    po::options_description count_opt("Error count options");
    count_opt.add_options()
    ("chrom-depth-file", po::value(&opt.chrom_depth_file),
     "If provided, the mean depth for each chromosome will be read from file, and these values will be used for high depth filtration. File should contain one line per chromosome, where each line begins with: \"chrom_name<TAB>depth\" (default: no chrom depth filtration)")
    ("max-depth-factor", po::value(&opt.max_depth_factor)->default_value(opt.max_depth_factor),
     "If a chrom depth file is supplied then loci with depth exceeding the mean chromosome depth times this value are filtered")
    ("counts-file", po::value(&opt.countsFilename),
     "write binary allele counts output to filename (required, no default)")
    ("observation-bed-file", po::value(&opt.observationsBedFilename),
     "write all observed indels to BED file (if not specified, individual indels will not be reported)")
    ("known-variants-vcf-file", po::value(&opt.knownVariantsFile),
     "VCF file specifying true variant genotypes.  Variants must be normalized.")
    ("excluded-regions-bed-file",
     po::value(&opt.excludedRegionsFileList)->multitoken(),
     "BED file specifying regions to exclude from counting. File must be tabix indexed. Argument can be provided multiple times to specify multiple exclusions.")
    ("nonempty-site-count-file",
     po::value(&opt.nonEmptySiteCountFilename),
     "File used to report the total number of non-empty sites observed which are otherwise eligible for error counting purposes. This file is used to monitor the approximate amount of evidence gathered")
    ;

    // final assembly
    po::options_description visible("Options");
    visible.add(aligndesc).add(count_opt);

    po::options_description visible2(get_starling_base_option_parser(opt));
    visible.add(visible2);

    po::options_description help_parse_opt("Help");
    help_parse_opt.add_options()
    ("help,h","print this message");

    visible.add(help_parse_opt);

    return visible;
}



void
finalizeSequenceAlleleCountsOptions(
    const prog_info& pinfo,
    const po::variables_map& vm,
    SequenceAlleleCountsOptions& opt)
{
    parseOptions(vm, opt.alignFileOpt);
    std::string errorMsg;
    if (checkOptions(opt.alignFileOpt, errorMsg))
    {
        pinfo.usage(errorMsg.c_str());
        //usage(log_os,prog,visible,errorMsg.c_str());
    }

    if (opt.countsFilename.empty())
    {
        pinfo.usage("Must specify a filename for allele counts output");
    }

    // knownVariantsFile and excludedRegionsFileList are both checked in the Python config code,
    // so we're not duplicating the effort here

    finalize_starling_base_options(pinfo,vm,opt);
}
