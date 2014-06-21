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

#include "blt_common/snp_pos_info.hh"
#include "blt_util/seq_util.hh"

#include <cassert>


/// \brief test if all basecalls in a pileup column are reference
///
inline
bool
is_spi_allref(const snp_pos_info& pi,
              const unsigned ref_gt)
{
    const unsigned n_calls(pi.calls.size());
    for (unsigned i(0); i<n_calls; ++i)
    {
        const uint8_t obs_id(pi.calls[i].base_id);
        assert(obs_id!=BASE_ID::ANY);
        if (ref_gt!=obs_id) return false;
    }
    return true;
}

