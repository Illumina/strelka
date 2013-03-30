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

/// \author Chris Saunders
///

#pragma once

#include "blt_common/snp_pos_info.hh"
#include "blt_util/seq_util.hh"

#include <cassert>


/// \brief test if all basecalls in a pileup column are reference
///
inline
bool
is_spi_allref(const snp_pos_info& pi,
              const unsigned ref_gt) {

    const unsigned n_calls(pi.calls.size());
    for(unsigned i(0); i<n_calls; ++i) {
        const uint8_t obs_id(pi.calls[i].base_id);
        assert(obs_id!=BASE_ID::ANY);
        if(ref_gt!=obs_id) return false;
    }
    return true;
}

