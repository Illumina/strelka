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

#include <cassert>

#include <algorithm>
#include <iosfwd>
#include <numeric>
#include <vector>


/// Discretizes support observations for a set of alleles at a single locus in a single sample
///
/// ultimately used for things like AD,ADF,ADR tags, or count-based EVS features
///
struct SupportingReadCountGroup
{
    SupportingReadCountGroup()
    {
        setAltCount(1);
    }

    /// set new alt count and clear
    void
    setAltCount(
        const unsigned altAlleleCount)
    {
        _confidentAlleleCount.resize(1+altAlleleCount);
        clear();
    }

    unsigned
    getAltCount() const
    {
        const unsigned s(_confidentAlleleCount.size());
        assert(s>0);
        return (s-1);
    }

    unsigned
    confidentRefAlleleCount() const
    {
        return _confidentAlleleCount[0];
    }

    unsigned
    confidentAltAlleleCount(const unsigned altAlleleIndex) const
    {
        return confidentAlleleCount(1+altAlleleIndex);
    }

    unsigned
    confidentAlleleCount(const unsigned alleleIndex) const
    {
        assert((alleleIndex) < _confidentAlleleCount.size());
        return _confidentAlleleCount[alleleIndex];
    }

    unsigned
    totalConfidentAltAlleleCount() const
    {
        return std::accumulate(_confidentAlleleCount.begin()+1,_confidentAlleleCount.end(),0);
    }

    unsigned
    totalConfidentCounts() const
    {
        return confidentRefAlleleCount() + totalConfidentAltAlleleCount();
    }

    void
    clear()
    {
        std::fill(_confidentAlleleCount.begin(),_confidentAlleleCount.end(),0);
        nonConfidentCount = 0;
    }

    void
    incrementRefAlleleCount()
    {
        incrementAlleleCount(0);
    }

    void
    incrementAltAlleleCount(const unsigned altAlleleIndex)
    {
        incrementAlleleCount(altAlleleIndex+1);
    }

    void
    incrementAlleleCount(
        const unsigned alleleIndex,
        const unsigned incrementBy = 1)
    {
        assert((alleleIndex) < _confidentAlleleCount.size());
        _confidentAlleleCount[alleleIndex] += incrementBy;
    }

    // number of ambiguous support reads
    unsigned nonConfidentCount = 0;
private:
    std::vector<unsigned> _confidentAlleleCount; ///< counts supporting alleles in REF,ALT1,ALT2.... order
};


/// Accumulate counts of confident supporting reads for each allele at a locus in one sample
///
/// this generalizes older indel support count data structures in that it is
/// designed with more than one alt allele in mind from the start
///
struct LocusSupportingReadStats
{
    void
    setAltCount(
        const unsigned altAlleleCount)
    {
        fwdCounts.setAltCount(altAlleleCount);
        revCounts.setAltCount(altAlleleCount);
    }

    unsigned
    getAltCount() const
    {
        return fwdCounts.getAltCount();
    }

    SupportingReadCountGroup&
    getCounts(const bool isFwdStrand)
    {
        return (isFwdStrand ? fwdCounts : revCounts);
    }

    const SupportingReadCountGroup&
    getCounts(const bool isFwdStrand) const
    {
        return (isFwdStrand ? fwdCounts : revCounts);
    }

    unsigned
    totalConfidentCounts() const
    {
        return getCounts(true).totalConfidentCounts() + getCounts(false).totalConfidentCounts();
    }

    void
    clear()
    {
        fwdCounts.clear();
        revCounts.clear();
    }

    SupportingReadCountGroup fwdCounts;
    SupportingReadCountGroup revCounts;
};


std::ostream&
operator<<(std::ostream& os, const LocusSupportingReadStats& lsrs);
