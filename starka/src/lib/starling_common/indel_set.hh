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

#ifndef __INDEL_SET_HH
#define __INDEL_SET_HH

#include "starling_common/indel.hh"

#include <iosfwd>
#include <set>


typedef std::set<indel_key> indel_set_t;


void
dump_indel_set(const indel_set_t& is,
               std::ostream& os);

#endif
