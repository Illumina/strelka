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
