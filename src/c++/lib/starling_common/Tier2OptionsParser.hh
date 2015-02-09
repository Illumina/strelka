// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
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

#include "Tier2Options.hh"

#include "boost/program_options.hpp"


boost::program_options::options_description
getTier2OptionsDescription(
    Tier2Options& opt);


bool
parseTier2Options(
    const boost::program_options::variables_map& vm,
    Tier2Options& opt,
    std::string& errorMsg);
