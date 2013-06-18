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

/// \author Chris Saunders
///

#ifndef __GVCF_HEADER
#define __GVCF_HEADER

#include "blt_common/blt_shared.hh"
#include "starling_common/chrom_depth_map.hh"

#include <iosfwd>
//#include <regex>


void
finish_gvcf_header(const gvcf_options& opt,
                   const cdmap_t& chrom_depth,
                   std::ostream& os);

std::string
determine_sample(const std::string bam_header_text);

#endif
