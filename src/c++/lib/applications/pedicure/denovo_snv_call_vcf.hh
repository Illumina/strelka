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
/// \author Morten Kallberg
///

#pragma once

#include "inovo_shared.hh"
#include "denovo_snv_call.hh"
#include "denovo_snv_call_info.hh"
#include "starling_common/PileupCleaner.hh"

#include <iosfwd>


void
denovo_snv_call_vcf(
    const inovo_options& opt,
    const inovo_deriv_options& dopt,
    const SampleInfoManager& sinfo,
    const cpiPtrTiers_t& pileups,
    const denovo_snv_call& dsc,
    std::ostream& os);
