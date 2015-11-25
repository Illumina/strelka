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

#include "snoise_option_parser.hh"

#include "../../starling_common/starling_base_option_parser.hh"
#include "blt_common/blt_arg_validate.hh"



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

    finalize_starling_base_options(pinfo,vm,opt);
}
