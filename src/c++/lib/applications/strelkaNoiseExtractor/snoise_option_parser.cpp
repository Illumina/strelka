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

#include "snoise_option_parser.hh"

#include "options/AlignmentFileOptionsParser.hh"
#include "starling_common/starling_base_option_parser.hh"



po::options_description
get_snoise_option_parser(
    snoise_options& opt)
{
    po::options_description aligndesc(getOptionsDescription(opt.alignFileOpt));

    po::options_description snoise_opt("Strelka noise extractor");
    snoise_opt.add_options()
    ("skip-vcf-header",
     po::value(&opt.is_skip_header)->zero_tokens(),
     "Skip vcf output header");

    po::options_description visible("Options");
    visible.add(aligndesc).add(snoise_opt);

    // add starling base options:
    po::options_description visible2(get_starling_base_option_parser(opt));
    visible.add(visible2);

    return visible;
}



void
finalize_snoise_options(
    const prog_info& pinfo,
    const po::variables_map& vm,
    snoise_options& opt)
{
    parseOptions(vm, opt.alignFileOpt);
    std::string errorMsg;
    if (checkOptions(opt.alignFileOpt, errorMsg))
    {
        pinfo.usage(errorMsg.c_str());
        //usage(log_os,prog,visible,errorMsg.c_str());
    }

    finalize_starling_base_options(pinfo,vm,opt);
}
