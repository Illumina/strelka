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

#pragma once


#include "blt_util/prog_info.hh"
#include "boost/program_options.hpp"
#include "starling_common/starling_base_shared.hh"

#include <iosfwd>


namespace po = boost::program_options;


/// test for the existence of an optional input file
///
/// filename has to either be empty or the file needs to exist and be accessible
///
/// if filename cannot be accessed this function will generate a usage error
///
/// label will be printed in any usage error message as "<label> file"
///
void
checkOptionalInputFile(
    const prog_info& pinfo,
    const std::string& filename,
    const char* label);


po::options_description
get_starling_base_option_parser(
    starling_base_options& opt);

/// validate options and process any required quick consistency
/// adjustments
///
void
finalize_starling_base_options(
    const prog_info& pinfo,
    const po::variables_map& vm,
    starling_base_options& opt);
