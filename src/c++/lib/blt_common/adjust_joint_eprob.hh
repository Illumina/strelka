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

/// \file
/// \author Chris Saunders
///

#pragma once

#include "blt_common/blt_shared.hh"
#include "blt_common/snp_pos_info.hh"

#include <vector>


/// cache dependent probs at a fixed vexp value (cache should be
/// associated with a specific blt_options object):
///
struct dependent_prob_cache
{
    enum { MAX_QSCORE = 64 };

    dependent_prob_cache() : _val(MAX_QSCORE+1), _is_init(MAX_QSCORE+1,false) {}

    // client is contracted to always call this function with the same
    // vexp value used on the first call
    //
    blt_float_t
    get_dependent_val(const unsigned qscore,
                      const blt_float_t vexp);

private:
    std::vector<blt_float_t> _val;
    std::vector<bool> _is_init;
};



void
adjust_joint_eprob(const blt_options& opt,
                   dependent_prob_cache& dpc,
                   const snp_pos_info& pi,
                   std::vector<float>& dependent_eprob);
