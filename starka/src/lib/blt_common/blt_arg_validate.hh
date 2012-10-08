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

#ifndef __BLT_ARG_VALIDATE_HH
#define __BLT_ARG_VALIDATE_HH


#include "blt_common/blt_shared.hh"
#include "blt_util/prog_info.hh"


extern const double MAX_MONOPLOID_THETA;
extern const double MAX_DIPLOID_THETA;


void
check_option_arg_range(const prog_info& pinfo,
                       const double val,
                       const char* label,
                       const double min,
                       const double max);


void
validate_blt_opt(const prog_info& pinfo,
                 const blt_options& client_opt,
                 const bool is_bacon_call_thresh,
                 const bool is_bacon_second_call_thresh,
                 const bool is_bacon_het_snp_ratio_thresh);

#endif
