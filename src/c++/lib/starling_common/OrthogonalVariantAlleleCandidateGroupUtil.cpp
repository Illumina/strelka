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

#include "OrthogonalVariantAlleleCandidateGroupUtil.hh"
#include "blt_util/prob_util.hh"

#include <vector>



template <typename T>
std::vector<size_t>
sortIndices(const std::vector<T>& v)
{
    // initialize original index locations
    std::vector<size_t> idx(v.size());
    for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

    // sort indexes based on comparing values in v
    std::sort(idx.begin(), idx.end(),
              [&v](size_t i1, size_t i2)
    {
        return v[i1] > v[i2];
    });

    return idx;
}



void
rankOrthogonalAllelesInSample(
    const unsigned sampleId,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    OrthogonalVariantAlleleCandidateGroup& rankedAlleleGroup,
    unsigned& referenceRank)
{
    const unsigned nonrefAlleleCount(alleleGroup.size());
    assert(nonrefAlleleCount!=0);

    // compile option below to use reads which only support a subset of alleles vs all alleles, no measurable difference
    // now, but could be important with longer alleles?
    std::set<unsigned> readIds;
#if 0
    // union of read ids
    for (unsigned nonrefAlleleIndex(0); nonrefAlleleIndex<nonrefAlleleCount; nonrefAlleleIndex++)
    {
        const IndelSampleData& isd(alleleGroup.data(nonrefAlleleIndex).getSampleData(sampleId));
        for (const auto& score : isd.read_path_lnp)
        {
            readIds.insert(score.first);
        }
    }
#else
    // intersection of read ids
    {
        std::map<unsigned,unsigned> countReadIds;
        for (unsigned nonrefAlleleIndex(0); nonrefAlleleIndex<nonrefAlleleCount; nonrefAlleleIndex++)
        {
            const IndelSampleData& isd(alleleGroup.data(nonrefAlleleIndex).getSampleData(sampleId));
            for (const auto& score : isd.read_path_lnp)
            {
                const auto iter(countReadIds.find(score.first));
                if (iter==countReadIds.end())
                {
                    countReadIds.insert(std::make_pair(score.first,1));
                }
                else
                {
                    iter->second += 1;
                }
            }
        }

        // filter countReadIds down to only the reads found for all alleles:
        for (const auto& value : countReadIds)
        {
            if (value.second >= nonrefAlleleCount)
            {
                readIds.insert(value.first);
            }
        }
    }
#endif

    // count of all haplotypes including reference
    const unsigned fullAlleleCount(nonrefAlleleCount+1);
    const unsigned refAlleleIndex(nonrefAlleleCount);

    std::vector<double> support(fullAlleleCount,0.);
    std::vector<double> lhood(fullAlleleCount);
    for (const auto readId : readIds)
    {
        bool isZeroAlleleCoverage(true);
        bool isPartialAlleleCoverage(false);
        for (unsigned nonrefAlleleIndex(0); nonrefAlleleIndex<nonrefAlleleCount; nonrefAlleleIndex++)
        {
            const IndelSampleData& isd(alleleGroup.data(nonrefAlleleIndex).getSampleData(sampleId));

            const auto iditer(isd.read_path_lnp.find(readId));
            if (iditer==isd.read_path_lnp.end())
            {
                isPartialAlleleCoverage=true;
                continue;
            }
            const ReadPathScores& path_lnp(iditer->second);

            if (isZeroAlleleCoverage)
            {
                lhood[refAlleleIndex] = static_cast<double>(path_lnp.ref);
            }
            else
            {
                lhood[refAlleleIndex] = std::max(lhood[refAlleleIndex],static_cast<double>(path_lnp.ref));
            }
            lhood[nonrefAlleleIndex] = path_lnp.indel;

            isZeroAlleleCoverage=false;
        }

        assert(not isZeroAlleleCoverage);

        // handle read which only supports a subset of alleles
        if (isPartialAlleleCoverage)
        {
            for (unsigned nonrefAlleleIndex(0); nonrefAlleleIndex < nonrefAlleleCount; nonrefAlleleIndex++)
            {
                const IndelSampleData& isd(alleleGroup.data(nonrefAlleleIndex).getSampleData(sampleId));

                const auto iditer(isd.read_path_lnp.find(readId));
                if (iditer != isd.read_path_lnp.end()) continue;

                lhood[nonrefAlleleIndex] = lhood[refAlleleIndex];
            }
        }

        unsigned maxIndex(0);
        normalize_ln_distro(lhood.begin(),lhood.end(),maxIndex);

        for (unsigned fullAlleleIndex(0); fullAlleleIndex<fullAlleleCount; fullAlleleIndex++)
        {
            support[fullAlleleIndex] += lhood[fullAlleleIndex];
        }
    }

    rankedAlleleGroup.clear();
    referenceRank = 0;

    bool isReferenceRankFound(false);
    for (const unsigned fullAlleleIndex : sortIndices(support))
    {
        if (fullAlleleIndex == refAlleleIndex)
        {
            isReferenceRankFound = true;
            continue;
        }
        if (not isReferenceRankFound) referenceRank++;
        rankedAlleleGroup.addVariantAllele(alleleGroup.iter(fullAlleleIndex));
    }
}



void
selectTopOrthogonalAllelesInSample(
    const unsigned sampleId,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const unsigned selectionSize,
    OrthogonalVariantAlleleCandidateGroup& topAlleleGroup)
{
    unsigned referenceRank(0);
    rankOrthogonalAllelesInSample(sampleId, alleleGroup, topAlleleGroup, referenceRank);

    unsigned topSize(selectionSize);
    if (referenceRank<topSize)
    {
        topSize -= 1;
    }

    if (topSize<alleleGroup.size())
    {
        topAlleleGroup.alleles.resize(topSize);
    }
}
