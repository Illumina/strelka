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

#include "inovo_shared.hh"
#include "starling_common/PileupCleaner.hh"

#include "blt_util/RegionProcessor.hh"


/// manage creation of a denovo callable bed track
struct DenovoCallableProcessor : public RegionProcessor
{
    typedef RegionProcessor base_t;

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
        using namespace INOVO_SAMPLETYPE;

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
