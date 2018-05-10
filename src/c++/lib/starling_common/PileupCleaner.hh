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

#include "blt_common/adjust_joint_eprob.hh"
#include "blt_common/blt_shared.hh"

#include <cassert>


/// filtered pileup with processed qualities and summary stats
struct CleanedPileup
{
    friend struct PileupCleaner;

    CleanedPileup()
        : _epi(cleanedPileup(),dependentErrorProb())
    {}

    unsigned
    totalBasecallCount() const
    {
        // pre-computed to reflect tier1/tier2
        return _n_raw_calls;
    }

    unsigned
    usedBasecallCount() const
    {
        return _cleanedPileup.calls.size();
    }

    unsigned
    unusedBasecallCount() const
    {
        return totalBasecallCount()- usedBasecallCount();
    }

    const snp_pos_info&
    rawPileup() const
    {
        assert(nullptr != _rawPileupPtr);
        return (*_rawPileupPtr);
    }

    const snp_pos_info&
    cleanedPileup() const
    {
        return _cleanedPileup;
    }

    const std::vector<float>&
    dependentErrorProb() const
    {
        return _dependentErrorProb;
    }

    /// deprecated, many old functions ask for this object, so this
    /// eases the transition
    const extended_pos_info&
    getExtendedPosInfo() const
    {
        return _epi;
    }

    void
    clear()
    {
        _rawPileupPtr = nullptr;
        _n_raw_calls = 0;
        _cleanedPileup.clear();
        _dependentErrorProb.clear();
    }

private:
    snp_pos_info&
    cleanedPileup()
    {
        return _cleanedPileup;
    }

    std::vector<float>&
    dependentErrorProb()
    {
        return _dependentErrorProb;
    }

    const snp_pos_info* _rawPileupPtr = nullptr;
    unsigned _n_raw_calls = 0;
    snp_pos_info _cleanedPileup;
    std::vector<float> _dependentErrorProb;
    const extended_pos_info _epi;
};


/// takes raw single sample pileup and processes information so that
/// it meets the criteria for snp calling
struct PileupCleaner
{
    explicit
    PileupCleaner(
        const blt_options& opt)
        : _opt(opt)
    {}

    void
    CleanPileupFilter(
        const snp_pos_info& pi,
        const bool is_include_tier2,
        CleanedPileup& cpi) const;

    void
    CleanPileupErrorProb(
        CleanedPileup& cpi) const;

    void
    CleanPileup(
        const snp_pos_info& pi,
        const bool is_include_tier2,
        CleanedPileup& cpi) const
    {
        CleanPileupFilter(pi, is_include_tier2, cpi);
        CleanPileupErrorProb(cpi);
    }

private:
    const blt_options& _opt;
    mutable dependent_prob_cache _dpcache;
};
