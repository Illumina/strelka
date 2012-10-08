// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
//

/// \file
///
/// \author Chris Saunders
///

#ifndef __STARLING_OPTION_PARSER_HH
#define __STARLING_OPTION_PARSER_HH

#include "blt_util/prog_info.hh"
#include "starling_common/starling_shared.hh"

#include "boost/program_options.hpp"

#include <iosfwd>


namespace po = boost::program_options;


// Called by program_options to parse avg_window_data:
//
void validate(boost::any& v,
              const std::vector<std::string>& values,
              std::vector<avg_window_data>*, int);


po::options_description
get_starling_shared_option_parser(starling_options& opt);

po::options_description
get_starling_option_parser(starling_options& opt);


void
write_starling_legacy_options(std::ostream& os);


// validate options and process any required quick consistency
// adjustments
//
void
finalize_starling_options(const prog_info& pinfo,
                          const po::variables_map& vm,
                          starling_options& opt);

#endif
