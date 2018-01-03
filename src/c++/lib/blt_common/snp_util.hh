//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
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

