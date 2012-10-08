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
