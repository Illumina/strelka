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

/// \file didn't have a better place to throw this function:

/// \author Chris Saunders
///

#include "blt_util/log.hh"
#include "blt_util/prob_util.hh"

#include <cstdlib>

#include <iomanip>
#include <iostream>



void
check_ln_distro_invalid_value(const char* label,
                              const double val,
                              const unsigned n)
{
    log_os << std::setprecision(14) << std::fixed;
    log_os << "ERROR: " << label << " element [" << n << "] has invalid value: '" << val << "'\n";
    log_os.unsetf(std::ios::fixed);
    exit(EXIT_FAILURE);
}



void
check_ln_distro_invalid_sum(const char* label,
                            const double sum)
{
    log_os << std::setprecision(14) << std::fixed;
    log_os << "ERROR: " << label << " sum is: '" << sum << "'\n";
    log_os.unsetf(std::ios::fixed);
    exit(EXIT_FAILURE);
}
