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

/// \author Chris Saunders
///

#pragma once

#include "inovo_shared.hh"
#include "denovo_indel_call.hh"
#include "starling_common/starling_indel_report_info.hh"

#include <array>
#include <iosfwd>


typedef std::array<starling_indel_sample_report_info,INOVO_TIERS::SIZE> isriTiers_t;


void
denovo_indel_call_vcf(
    const inovo_options& opt,
    const inovo_deriv_options& dopt,
    const SampleInfoManager& sinfo,
    const denovo_indel_call& dinc,
    const starling_indel_report_info& iri,
    const std::vector<isriTiers_t>& isri,
    std::ostream& os);
