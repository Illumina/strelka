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
#include "indel_util.hh"
#include "blt_util/prob_util.hh"

#include <vector>

// turn this on to use reads which support some, but not all, of an
// overlapping allele group;
// #define USE_GERMLINE_SUPPORTING_READ_UNION



void
getAlleleGroupUnionReadIds(
    const unsigned sampleId,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    std::set<unsigned>& readIds,
    const bool is_tier1_only)
{
    const unsigned nonrefAlleleCount(alleleGroup.size());
    for (unsigned nonrefAlleleIndex(0); nonrefAlleleIndex < nonrefAlleleCount; nonrefAlleleIndex++)
    {
        const IndelSampleData& isd(alleleGroup.data(nonrefAlleleIndex).getSampleData(sampleId));
        for (const auto& score : isd.read_path_lnp)
        {
            if (is_tier1_only && (! score.second.is_tier1_read)) continue;

            readIds.insert(score.first);
        }
    }
}



void
getAlleleGroupIntersectionReadIds(
    const unsigned sampleId,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    std::set<unsigned>& readIds,
    const bool isTier1Only)
{
    const unsigned nonrefAlleleCount(alleleGroup.size());
    std::map<unsigned,unsigned> countReadIds;
    for (unsigned nonrefAlleleIndex(0); nonrefAlleleIndex<nonrefAlleleCount; nonrefAlleleIndex++)
    {
        const IndelSampleData& isd(alleleGroup.data(nonrefAlleleIndex).getSampleData(sampleId));
        for (const auto& score : isd.read_path_lnp)
        {
            if (isTier1Only && (! score.second.is_tier1_read)) continue;

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



/// find set of read ids which support the entire set of alleles in alleleGroup
///
static
void
getAlleleGroupSupportingReadIds(
    const unsigned sampleId,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    std::set<unsigned>& readIds)
{
    static const bool isTier1Only(false);

#ifdef USE_GERMLINE_SUPPORTING_READ_UNION
    getAlleleGroupUnionReadIds(sampleId, alleleGroup, readIds, isTier1Only);
#else
    getAlleleGroupIntersectionReadIds(sampleId, alleleGroup, readIds, isTier1Only);
#endif
}



void
getAlleleLikelihoodsFromRead(
    const unsigned sampleId,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const unsigned readId,
    std::vector<double>& lhood)
{
    const unsigned nonrefAlleleCount(alleleGroup.size());
    const unsigned fullAlleleCount(nonrefAlleleCount+1);
    const unsigned refAlleleIndex(nonrefAlleleCount);

    lhood.resize(fullAlleleCount);

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
}



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

    std::set<unsigned> readIds;
    getAlleleGroupSupportingReadIds(sampleId, alleleGroup, readIds);

    // count of all haplotypes including reference
    const unsigned fullAlleleCount(nonrefAlleleCount+1);
    const unsigned refAlleleIndex(nonrefAlleleCount);

    std::vector<double> support(fullAlleleCount,0.);
    for (const auto readId : readIds)
    {
        std::vector<double> lhood(fullAlleleCount);
        getAlleleLikelihoodsFromRead(sampleId, alleleGroup, readId, lhood);
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



bool
addAllelesAtOtherPositions(
    const pos_t pos,
    const pos_t largest_total_indel_ref_span_per_read,
    const unsigned sampleIndex,
    const int callerPloidy,
    const IndelBuffer& indelBuffer,
    OrthogonalVariantAlleleCandidateGroup& alleleGroup)
{
    const pos_t minIndelBufferPos(pos-largest_total_indel_ref_span_per_read);

    const unsigned inputAlleleCount(alleleGroup.size());
    assert(inputAlleleCount!=0);

    bool isEveryAltOrthogonal(true);

    // first get the set of candidate alt alleles from positions other than 'pos'
    //
    std::vector<IndelKey> filteredAltAlleles;
    {
        const known_pos_range inputAlleleGroupRange(alleleGroup.getReferenceRange());

        // extend end_pos by one to ensure that we find indels adjacent to the right end of the range
        /// TODO add strict definition and unit tests to rangeIterator wrt adjacent indels
        const auto indelIterPair(indelBuffer.rangeIterator(inputAlleleGroupRange.begin_pos, inputAlleleGroupRange.end_pos+1));
        for (auto altAlleleIter(indelIterPair.first); altAlleleIter!=indelIterPair.second; ++altAlleleIter)
        {
            const IndelKey& altAlleleKey(altAlleleIter->first);

            // all alleles with this starting position have already been
            // considered..
            if (altAlleleKey.pos == pos) continue;

            // filter out indels which have already been cleared out of the indel buffer:
            if (altAlleleKey.pos < minIndelBufferPos) continue;

            // no breakpoints:
            if (altAlleleKey.is_breakpoint()) continue;

            // must be a candidate allele:
            const IndelData& altAlleleData(getIndelData(altAlleleIter));
            if (not indelBuffer.isCandidateIndel(altAlleleKey, altAlleleData)) continue;

            // must be orthogonal to all input alleles:
            {
                bool isOrthogonalToAllInputAlleles(true);
                for (unsigned inputAlleleIndex(0); inputAlleleIndex < inputAlleleCount; inputAlleleIndex++)
                {
                    const IndelKey& inputAlleleKey(alleleGroup.key(inputAlleleIndex));
                    if (is_indel_conflict(altAlleleKey, inputAlleleKey)) continue;
                    isOrthogonalToAllInputAlleles = false;
                    break;
                }
                if (not isOrthogonalToAllInputAlleles)
                {
                    isEveryAltOrthogonal = false;
                    continue;
                }
            }

            // made it! add allele to the set we move forward with:
            filteredAltAlleles.push_back(altAlleleKey);
        }
    }

    if (filteredAltAlleles.empty()) return isEveryAltOrthogonal;

    OrthogonalVariantAlleleCandidateGroup altAlleleGroup;
    for (const auto& altAlleleKey : filteredAltAlleles)
    {
        const auto altIter(indelBuffer.getIndelIter(altAlleleKey));
        altAlleleGroup.addVariantAllele(altIter);
    }

    if (altAlleleGroup.size()>1)
    {
        // rank alt alleles and include from highest to lowest unless interference clique is broken:
        unsigned referenceRank(0);
        OrthogonalVariantAlleleCandidateGroup rankedAltAlleleGroup;
        rankOrthogonalAllelesInSample(sampleIndex, altAlleleGroup, rankedAltAlleleGroup, referenceRank);

        altAlleleGroup.clear();

        for (const auto& rankedAltAlleleIter : rankedAltAlleleGroup.alleles)
        {
            bool isGroupOrthogonal(true);
            for (const auto& altAlleleIter : altAlleleGroup.alleles)
            {
                if (not is_indel_conflict(rankedAltAlleleIter->first, altAlleleIter->first))
                {
                    isGroupOrthogonal = false;
                    break;
                }
            }
            if (isGroupOrthogonal)
            {
                altAlleleGroup.alleles.push_back(rankedAltAlleleIter);
            }
        }
    }

    // put all qualifying alts back together with variants to form an extended allele set:
    OrthogonalVariantAlleleCandidateGroup extendedVariantAlleleGroup(alleleGroup);
    for (const auto& altAlleleIter : altAlleleGroup.alleles)
    {
        extendedVariantAlleleGroup.alleles.push_back(altAlleleIter);
    }

    // rerank and reselect top N alleles, N=callerPloidy
    //
    selectTopOrthogonalAllelesInSample(sampleIndex, extendedVariantAlleleGroup, callerPloidy, alleleGroup);

    return isEveryAltOrthogonal;
}
