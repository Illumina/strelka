// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

/// \file
///
/// \author Chris Saunders
///

#ifndef __STRELKA_OPTION_PARSER_HH
#define __STRELKA_OPTION_PARSER_HH

#include "blt_util/prog_info.hh"
#include "strelka/strelka_shared.hh"

#include "boost/program_options.hpp"


namespace po = boost::program_options;


po::options_description
get_strelka_option_parser(strelka_options& opt);


// validate options and process any required quick consistency
// adjustments
//
void
finalize_strelka_options(const prog_info& pinfo,
                         const po::variables_map& vm,
                         strelka_options& opt);

#endif
