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

#include "PileupCleaner.hh"



void
PileupCleaner::
CleanPileupFilter(
    const snp_pos_info& pi,
    const bool is_include_tier2,
    CleanedPileup& cpi) const
{
    cpi.clear();
    cpi._rawPileupPtr = &pi;
    snp_pos_info& cleanedPi(cpi.cleanedPileup());

    cleanedPi.set_ref_base(pi.get_ref_base());

    cpi._n_raw_calls = pi.calls.size();
    for (const auto& bc : pi.calls)
    {
        if (bc.is_call_filter)
        {
            if (! (is_include_tier2 &&
                   bc.is_tier_specific_call_filter))
            {
                continue;
            }
        }
        cleanedPi.calls.push_back(bc);
    }

    if (is_include_tier2)
    {
        cpi._n_raw_calls += pi.tier2_calls.size();
        for (const auto& bc : pi.tier2_calls)
        {
            if (bc.is_call_filter) continue;
            cleanedPi.calls.push_back(bc);
        }
    }
}



void
PileupCleaner::
CleanPileupErrorProb(
    CleanedPileup& cpi) const
{
    adjust_joint_eprob(_opt, _dpcache, cpi.cleanedPileup(), cpi.dependentErrorProb());
}
