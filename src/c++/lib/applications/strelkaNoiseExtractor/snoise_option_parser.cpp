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

#include "snoise_option_parser.hh"
#include "blt_common/blt_arg_validate.hh"
#include "starling_common/starling_option_parser.hh"



po::options_description
get_snoise_option_parser(
    snoise_options& opt)
{
    po::options_description snoise_opt("Strelka noise extractor");
    snoise_opt.add_options()
    ("skip-vcf-header",
     po::value(&opt.is_skip_header)->zero_tokens(),
     "Skip vcf output header");

    po::options_description visible("Options");
    visible.add(snoise_opt);

    return visible;
}



void
finalize_snoise_options(
    const prog_info& pinfo,
    const po::variables_map& vm,
    snoise_options& opt)
{
    // nothing to validate in this layer:

    finalize_starling_options(pinfo,vm,opt);
}
