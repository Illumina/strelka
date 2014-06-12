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
