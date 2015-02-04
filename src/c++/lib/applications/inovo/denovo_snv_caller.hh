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

#include "denovo_snv_call.hh"
#include "inovo_shared.hh"
#include "blt_util/prog_info.hh"
#include "starling_common/PileupCleaner.hh"

void
get_denovo_snv_call(
    const inovo_options& opt,
    const SampleInfoManager& sinfo,
    const std::vector<const CleanedPileup*>& pileups,
    denovo_snv_call& dsc);
