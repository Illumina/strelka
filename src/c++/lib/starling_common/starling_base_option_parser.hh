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

#pragma once


#include "blt_util/prog_info.hh"
#include "boost/program_options.hpp"
#include "starling_common/starling_base_shared.hh"

#include <iosfwd>


namespace po = boost::program_options;


// Called by program_options to parse avg_window_data:
//
void validate(
    boost::any& v,
    const std::vector<std::string>& values,
    std::vector<avg_window_data>*, int);


po::options_description
get_starling_base_option_parser(
    starling_base_options& opt);

void
write_starling_legacy_options(
    const starling_base_options& default_opt,
    std::ostream& os);


/// validate options and process any required quick consistency
/// adjustments
///
void
finalize_starling_base_options(
    const prog_info& pinfo,
    const po::variables_map& vm,
    starling_base_options& opt);
