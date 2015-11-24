// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
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

#include "pedicure_shared.hh"
#include "starling_common/PileupCleaner.hh"

#include "blt_util/RegionProcessor.hh"


/// manage creation of a denovo callable bed track
struct DenovoCallableProcessor : public RegionProcessor
{
    typedef RegionProcessor base_t;

    explicit
    DenovoCallableProcessor(
        std::ostream* osptr) :
        base_t(osptr)
    {}

    void
    addToRegion(
        const std::string& chrom,
        const pos_t outputPos,
        const SampleInfoManager& sinfo,
        const std::vector<const CleanedPileup*>& pileups)
    {
        using namespace PEDICURE_SAMPLETYPE;

        const unsigned probandIndex(sinfo.getTypeIndexList(PROBAND)[0]);
        if (testIndex(pileups,probandIndex)) return;

        const std::vector<unsigned>& parentIndex(sinfo.getTypeIndexList(PARENT));
        for (const auto sampleIndex : parentIndex)
        {
            if (testIndex(pileups,sampleIndex)) return;
        }

        base_t::addToRegion(chrom,outputPos);
    }

private:
    bool
    testIndex(
        const std::vector<const CleanedPileup*>& pileups,
        const unsigned index) const
    {
        static const unsigned minDepth(18);

        const CleanedPileup& cpi(*pileups[index]);
        return (cpi.cleanedPileup().calls.size() < minDepth);
    }
};
