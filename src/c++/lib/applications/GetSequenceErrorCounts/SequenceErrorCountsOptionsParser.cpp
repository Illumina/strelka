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

#include "SequenceErrorCountsOptionsParser.hh"

#include "boost/filesystem.hpp"



//#define DEBUG_OPTIONS
#ifdef DEBUG_OPTIONS
#include "blt_util/log.hh"
#endif




po::options_description
getSequenceErrorCountsOptionsParser(
    SequenceErrorCountsOptions& opt)
{
    po::options_description count_opt("Error count options");
    count_opt.add_options()
    ("chrom-depth-file", po::value(&opt.chrom_depth_file),
     "If provided, the mean depth for each chromosome will be read from file, and these values will be used for high depth filtration. File should contain one line per chromosome, where each line begins with: \"chrom_name<TAB>depth\" (default: no chrom depth filtration)")
    ("max-depth-factor", po::value(&opt.max_depth_factor)->default_value(opt.max_depth_factor),
     "If a chrom depth file is supplied then loci with depth exceeding the mean chromosome depth times this value are filtered")
    ("counts-file", po::value(&opt.countsFilename),
     "write binary error counts output to filename (required, no default)")
    ("observation-bed-file", po::value(&opt.observationsBedFilename),
     "write all observed indels to BED file (if not specified, individual indels will not be reported)")
    ;

    // final assembly
    po::options_description visible("Options");
    visible.add(count_opt);

    po::options_description visible2(get_starling_base_option_parser(opt));
    visible.add(visible2);

    po::options_description help_parse_opt("Help");
    help_parse_opt.add_options()
    ("help,h","print this message");

    visible.add(help_parse_opt);

    return visible;
}



void
finalizeSequenceErrorCountsOptions(
    const prog_info& pinfo,
    const po::variables_map& vm,
    SequenceErrorCountsOptions& opt)
{
    if (opt.bam_filename.empty())
    {
        pinfo.usage("Must specify a sorted & indexed BAM/CRAM file containing aligned sample reads");
    }
    else
    {
        if (! boost::filesystem::exists(opt.bam_filename))
        {
            std::ostringstream oss;
            oss << "Submitted BAM/CRAM file does not exist: '" << opt.bam_filename << "'";
            pinfo.usage(oss.str().c_str());
        }
    }

    if (opt.countsFilename.empty())
    {
        pinfo.usage("Must specify an filename for error counts output");
    }


    /// not using any of these below options, not worth cleaning this up right now...

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

    if (! opt.germline_variant_scoring_models_filename.empty())
    {
        if (! boost::filesystem::exists(opt.germline_variant_scoring_models_filename))
        {
            std::ostringstream oss;
            oss << "Germline scoring model file does not exist: '" << opt.germline_variant_scoring_models_filename << "'";
            pinfo.usage(oss.str().c_str());
        }
    }
    else
    {
        if (! opt.germline_variant_scoring_model_name.empty())
        {
            pinfo.usage("Specified germline variant scoring model name without a corresponding scoring file.");
        }
    }

    finalize_starling_base_options(pinfo,vm,opt);
}

