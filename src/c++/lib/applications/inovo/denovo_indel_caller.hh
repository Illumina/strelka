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

#include "denovo_indel_call.hh"
#include "inovo_shared.hh"


void
get_denovo_indel_call(
    const inovo_options& opt,
    const SampleInfoManager& sinfo,
    const std::vector<const starling_sample_options*>& sampleOptions,
    const double indel_error_prob,
    const double ref_error_prob,
    const indel_key& ik,
    const std::vector<const indel_data*>& allIndelData,
    const bool is_use_alt_indel,
    denovo_indel_call& dinc);
