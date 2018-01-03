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

///
/// \author Chris Saunders
///

#pragma once

#include "strelka_shared.hh"

#include "blt_util/prog_info.hh"

#include "boost/program_options.hpp"


namespace po = boost::program_options;


po::options_description
get_strelka_option_parser(
    strelka_options& opt);


// validate options and process any required quick consistency
// adjustments
//
void
finalize_strelka_options(
    const prog_info& pinfo,
    const po::variables_map& vm,
    strelka_options& opt);
