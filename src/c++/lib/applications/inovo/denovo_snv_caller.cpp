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

/// \author Chris Saunders
///

#include "denovo_snv_caller.hh"

#include <array>



static
bool
checkME(
    const uint8_t a1,
    const uint8_t a2,
    const std::array<unsigned,N_BASE>& p1,
    const std::array<unsigned,N_BASE>& p2)
{
    static const unsigned minCount(2);
    return ((p1[a1]>=minCount) && (p2[a2]>=minCount));
}



void
get_denovo_snv_call(
    const inovo_options& opt,
    const SampleInfoManager& sinfo,
    const std::vector<const CleanedPileup*>& pileups,
    denovo_snv_call& dsc)
{
    // simple count based rules:
    // depth must be at least 10 in all of proband and parents
    // enumerate alleles in proband, take the top two which occur at a rate higher than 20%
    // if two alleles, find at least 2 observations of each allele coming from each parent
    // if one allele, find at least 2 observation of the allele coming from each parent
    //

    using namespace INOVO_SAMPLETYPE;

    const unsigned probandIndex(sinfo.getTypeIndexList(PROBAND)[0]);
    const std::vector<unsigned>& parentIndex(sinfo.getTypeIndexList(PARENT));

    std::vector<unsigned> allIndex(parentIndex);
    allIndex.push_back(probandIndex);

    for (const auto sampleIndex : allIndex)
    {
        const CleanedPileup& cpi(*pileups[sampleIndex]);
        if (cpi.cleanedPileup().calls.size() < 10)
        {
            return;
        }
    }

    // top two alleles;
    uint8_t max1(0), max2(0);
    {
        const CleanedPileup& probandCpi(*pileups[probandIndex]);

        std::array<unsigned,N_BASE> base_count;
        probandCpi.cleanedPileup().get_known_counts(base_count,opt.min_qscore);

        unsigned total(0);
        for (const unsigned count : base_count)
        {
            total += count;
        }
        const unsigned min_allele_count(total*.2);

        for (unsigned i(0);i<N_BASE;++i)
        {
            if(base_count[i] > base_count[max1]) max1=i;
        }

        bool is_max2(false);
        for (unsigned i(0);i<N_BASE;++i)
        {
            if (i==max1) continue;
            if (is_max2 && (base_count[i] <= base_count[max2])) continue;
            if (base_count[i]<min_allele_count) continue;
            max2=i;
            is_max2=true;
        }

        if (! is_max2) max2 = max1;
    }


    // the nature of the hack count model requires two parents:
    if (parentIndex.size() < 2) return;

    //
    {
        const CleanedPileup& p1Cpi(*pileups[parentIndex[0]]);
        const CleanedPileup& p2Cpi(*pileups[parentIndex[1]]);

        std::array<unsigned,N_BASE> p1Count;
        p1Cpi.cleanedPileup().get_known_counts(p1Count,opt.min_qscore);

        std::array<unsigned,N_BASE> p2Count;
        p2Cpi.cleanedPileup().get_known_counts(p2Count,opt.min_qscore);

        // scenario 1: max1 -> parent1, max2 -> parent2
        bool isValid(checkME(max1,max2,p1Count,p2Count));
        if (! isValid)
        {
            // scenario 2: max1 -> parent2, max2 -> parent1
            isValid = checkME(max1,max2,p2Count,p1Count);
        }
        if (! isValid)
        {
            return;
        }
    }

    dsc.rs.isCall = true;
}
