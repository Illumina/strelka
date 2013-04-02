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
                 const blt_options& client_opt);

#endif
