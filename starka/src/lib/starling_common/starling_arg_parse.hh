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
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the begining of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///

#ifndef __STARLING_ARG_PARSE_H
#define __STARLING_ARG_PARSE_H

#include "blt_common/blt_arg_parse_util.hh"
#include "starling_common/starling_shared.hh"


// this combines starling argument parsing and validation:
//
void
legacy_starling_arg_parse(arg_data& ad,
                          starling_options& client_opt);

#endif
