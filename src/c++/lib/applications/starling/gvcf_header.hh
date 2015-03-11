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

///
/// \author Chris Saunders
///

#pragma once

#include "calibration_models.hh"
#include "starling_shared.hh"
#include "blt_common/blt_shared.hh"
#include "blt_util/chrom_depth_map.hh"

#include <iosfwd>


void
finish_gvcf_header(const starling_options& opt,
                   const gvcf_deriv_options& dopt,
                   const cdmap_t& chrom_depth,
                   const std::string& bam_header_data,
                   std::ostream& os,
                   calibration_models& CM);
