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

#include "starling_shared.hh"
#include "starling_common/starling_base_option_parser.hh"


po::options_description
get_starling_option_parser(starling_options& opt);


// validate options and process any required quick consistency
// adjustments
//
void
finalize_starling_options(
    const prog_info& pinfo,
    const po::variables_map& vm,
    starling_options& opt);

