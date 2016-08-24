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
#include "blt_util/sort_util.hh"

#include <vector>

// turn this on to use reads which support some, but not all, of an
// overlapping allele group;
// #define USE_GERMLINE_SUPPORTING_READ_UNION



void
getAlleleGroupUnionReadIds(
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    std::set<unsigned>& readIds,
    const bool isTier1Only)
{
    const unsigned nonrefAlleleCount(alleleGroup.size());
    for (unsigned nonrefAlleleIndex(0); nonrefAlleleIndex < nonrefAlleleCount; nonrefAlleleIndex++)
    {
        const IndelSampleData& isd(alleleGroup.data(nonrefAlleleIndex).getSampleData(sampleIndex));
        for (const auto& score : isd.read_path_lnp)
        {
            if (isTier1Only && (! score.second.is_tier1_read)) continue;

            readIds.insert(score.first);
        }
    }
}



void
getAlleleGroupIntersectionReadIds(
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    std::set<unsigned>& readIds,
    const bool isTier1Only)
{
    const unsigned nonrefAlleleCount(alleleGroup.size());
    std::map<unsigned,unsigned> countReadIds;
    for (unsigned nonrefAlleleIndex(0); nonrefAlleleIndex<nonrefAlleleCount; nonrefAlleleIndex++)
    {
        const IndelSampleData& isd(alleleGroup.data(nonrefAlleleIndex).getSampleData(sampleIndex));
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



void
getAlleleGroupSupportingReadIds(
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    std::set<unsigned>& readIds,
    const bool isTier1Only)
{
#ifdef USE_GERMLINE_SUPPORTING_READ_UNION
    getAlleleGroupUnionReadIds(sampleId, alleleGroup, readIds, isTier1Only);
#else
    getAlleleGroupIntersectionReadIds(sampleIndex, alleleGroup, readIds, isTier1Only);
#endif
}



void
getAlleleLogLhoodFromRead(
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const unsigned readId,
    std::vector<double>& alleleLogLhood)
{
    const unsigned nonrefAlleleCount(alleleGroup.size());
    const unsigned fullAlleleCount(nonrefAlleleCount+1);
    const unsigned refAlleleIndex(nonrefAlleleCount);

    alleleLogLhood.resize(fullAlleleCount);

    bool isZeroAlleleCoverage(true);
    bool isPartialAlleleCoverage(false);
    for (unsigned nonrefAlleleIndex(0); nonrefAlleleIndex<nonrefAlleleCount; nonrefAlleleIndex++)
    {
        const IndelSampleData& isd(alleleGroup.data(nonrefAlleleIndex).getSampleData(sampleIndex));

        const auto iditer(isd.read_path_lnp.find(readId));
        if (iditer==isd.read_path_lnp.end())
        {
            isPartialAlleleCoverage=true;
            continue;
        }
        const ReadPathScores& path_lnp(iditer->second);

        if (isZeroAlleleCoverage)
        {
            alleleLogLhood[refAlleleIndex] = static_cast<double>(path_lnp.ref);
        }
        else
        {
            alleleLogLhood[refAlleleIndex] = std::max(alleleLogLhood[refAlleleIndex],static_cast<double>(path_lnp.ref));
        }
        alleleLogLhood[nonrefAlleleIndex] = path_lnp.indel;

        isZeroAlleleCoverage=false;
    }

    assert(not isZeroAlleleCoverage);

    // handle read which only supports a subset of alleles
    if (isPartialAlleleCoverage)
    {
        for (unsigned nonrefAlleleIndex(0); nonrefAlleleIndex < nonrefAlleleCount; nonrefAlleleIndex++)
        {
            const IndelSampleData& isd(alleleGroup.data(nonrefAlleleIndex).getSampleData(sampleIndex));

            const auto iditer(isd.read_path_lnp.find(readId));
            if (iditer != isd.read_path_lnp.end()) continue;

            alleleLogLhood[nonrefAlleleIndex] = alleleLogLhood[refAlleleIndex];
        }
    }
}



void
getAlleleNaivePosteriorFromRead(
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const unsigned readId,
    std::vector<double>& alleleProb)
{
    getAlleleLogLhoodFromRead(sampleIndex, alleleGroup, readId, alleleProb);
    unsigned maxIndex(0);
    normalize_ln_distro(alleleProb.begin(),alleleProb.end(),maxIndex);
}


void
rankOrthogonalAllelesInSample(
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    OrthogonalVariantAlleleCandidateGroup& rankedAlleleGroup,
    unsigned& referenceRank)
{
    const unsigned nonrefAlleleCount(alleleGroup.size());
    assert(nonrefAlleleCount!=0);

    static const bool isTier1Only(false);
    std::set<unsigned> readIds;
    getAlleleGroupSupportingReadIds(sampleIndex, alleleGroup, readIds, isTier1Only);

    // count of all haplotypes including reference
    const unsigned fullAlleleCount(nonrefAlleleCount+1);
    const unsigned refAlleleIndex(nonrefAlleleCount);

    std::vector<double> support(fullAlleleCount,0.);
    for (const auto readId : readIds)
    {
        std::vector<double> lhood(fullAlleleCount);
        getAlleleNaivePosteriorFromRead(sampleIndex, alleleGroup, readId, lhood);
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
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& inputAlleleGroup,
    const unsigned selectionSize,
    OrthogonalVariantAlleleCandidateGroup& topAlleleGroup)
{
    unsigned referenceRank(0);
    rankOrthogonalAllelesInSample(sampleIndex, inputAlleleGroup, topAlleleGroup, referenceRank);

    unsigned topSize(selectionSize);
    if (referenceRank<topSize)
    {
        topSize -= 1;
    }

    if (topSize<inputAlleleGroup.size())
    {
        topAlleleGroup.alleles.resize(topSize);
    }
}



void
selectTopOrthogonalAllelesInAllSamples(
    const unsigned sampleCount,
    const std::vector<unsigned>& callerPloidy,
    const OrthogonalVariantAlleleCandidateGroup& inputAlleleGroup,
    OrthogonalVariantAlleleCandidateGroup& topAlleleGroup,
    std::vector<unsigned>& topVariantAlleleIndexPerSample)
{
    assert(sampleCount == callerPloidy.size());

    // this structure is used to create an approximate allele rank over all samples
    std::map<IndelKey,unsigned> topVariantAlleleKeyScore;

    std::vector<IndelKey> topVariantAllelePerSample;

    for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
    {
        // in each sample, rank and select top N alleles, N=callerPloidy, and accumulate these top alleles
        // over all samples in topAlleleGroup
        //
        const unsigned sampleCallerPloidy(callerPloidy[sampleIndex]);
        OrthogonalVariantAlleleCandidateGroup topAlleleGroupInSample;
        assert(callerPloidy[sampleIndex] > 0);
        selectTopOrthogonalAllelesInSample(sampleIndex, inputAlleleGroup, sampleCallerPloidy,
                                           topAlleleGroupInSample);

        const unsigned topAllelesInSampleCount(topAlleleGroupInSample.size());
        for (unsigned alleleIndex(0); alleleIndex < topAllelesInSampleCount; alleleIndex++)
        {
            const IndelKey& indelKey(topAlleleGroupInSample.key(alleleIndex));
            auto indelKeyIter(topVariantAlleleKeyScore.find(indelKey));
            if (indelKeyIter == topVariantAlleleKeyScore.end())
            {
                auto retVal = topVariantAlleleKeyScore.insert(std::make_pair(indelKey,0));
                indelKeyIter = retVal.first;
                topAlleleGroup.addVariantAllele(topAlleleGroupInSample.iter(alleleIndex));
            }
            indelKeyIter->second += (sampleCallerPloidy-alleleIndex);

            if (alleleIndex == 0)
            {
                topVariantAllelePerSample[sampleIndex] = indelKey;
            }
        }
    }

    // approximately rank topVariantAlleleGroup alleles based on sample rankings
    if (sampleCount > 1)
    {
        std::vector<int> support;
        const unsigned topAllelesCount(topAlleleGroup.size());
        for (unsigned alleleIndex(0); alleleIndex < topAllelesCount; alleleIndex++)
        {
            const IndelKey& indelKey(topAlleleGroup.key(alleleIndex));
            const auto indelKeyIter(topVariantAlleleKeyScore.find(indelKey));
            assert(indelKeyIter != topVariantAlleleKeyScore.end());
            support.push_back(indelKeyIter->second);
        }

        const OrthogonalVariantAlleleCandidateGroup tmpCopy(topAlleleGroup);
        topAlleleGroup.clear();
        for (const unsigned altAlleleIndex : sortIndices(support))
        {
            topAlleleGroup.addVariantAllele(tmpCopy.iter(altAlleleIndex));
        }
    }

    // convert top variant per sample from indel key to topAlleleGroup index:
    {
        std::map<IndelKey,unsigned> indelIndex;
        const unsigned topAllelesCount(topAlleleGroup.size());
        for (unsigned alleleIndex(0); alleleIndex < topAllelesCount; alleleIndex++)
        {
            const IndelKey& indelKey(topAlleleGroup.key(alleleIndex));
            indelIndex[indelKey] = alleleIndex;
        }
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            const auto indelIter(indelIndex.find(topVariantAllelePerSample[sampleIndex]));
            assert(indelIter != indelIndex.end());
            topVariantAlleleIndexPerSample[sampleIndex] = indelIter->second;
        }
    }
}



bool
addAllelesAtOtherPositions(
    const unsigned sampleCount,
    const std::vector<unsigned>& callerPloidy,
    const pos_t pos,
    const pos_t largest_total_indel_ref_span_per_read,
    const IndelBuffer& indelBuffer,
    OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    std::vector<unsigned>& topVariantAlleleIndexPerSample)
{
    const pos_t minIndelBufferPos(pos-largest_total_indel_ref_span_per_read);

    const unsigned inputAlleleCount(alleleGroup.size());
    assert(inputAlleleCount!=0);

    bool isEveryAltOrthogonal(true);

    // first get the set of candidate alt alleles from positions other than 'pos'
    //
    OrthogonalVariantAlleleCandidateGroup newAltAlleleGroup;
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
            const auto altIter(indelBuffer.getIndelIter(altAlleleKey));
            newAltAlleleGroup.addVariantAllele(altIter);
        }
    }

    if (newAltAlleleGroup.empty()) return isEveryAltOrthogonal;


    const unsigned newAltAlleleCount(newAltAlleleGroup.size());
    if (newAltAlleleCount > 1)
    {
        // rank new alt alleles according to support level over all samples
        OrthogonalVariantAlleleCandidateGroup rankedNewAltAlleleGroup;
        {
            // this structure is used to create an approximate aggregate allele rank over all samples
            std::map<IndelKey, unsigned> newAltAlleleKeyScore;
            for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
            {
                // rank alt alleles and include from highest to lowest unless interference clique is broken:
                //
                /// TODO shouldn't this ranking be done wrt all alleles and not just the new alts?
                unsigned referenceRank(0);
                OrthogonalVariantAlleleCandidateGroup rankedNewAltAlleleGroupPerSample;
                rankOrthogonalAllelesInSample(sampleIndex, newAltAlleleGroup, rankedNewAltAlleleGroupPerSample,
                                              referenceRank);

                unsigned refPenalty(0);
                for (unsigned alleleIndex(0); alleleIndex < newAltAlleleCount; alleleIndex++)
                {
                    const IndelKey& indelKey(rankedNewAltAlleleGroupPerSample.key(alleleIndex));
                    auto indelKeyIter(newAltAlleleKeyScore.find(indelKey));
                    if (indelKeyIter == newAltAlleleKeyScore.end())
                    {
                        auto retVal = newAltAlleleKeyScore.insert(std::make_pair(indelKey, 0));
                        indelKeyIter = retVal.first;
                        rankedNewAltAlleleGroup.addVariantAllele(rankedNewAltAlleleGroupPerSample.iter(alleleIndex));

                    }
                    if (referenceRank == alleleIndex) refPenalty = 1;
                    indelKeyIter->second += ((newAltAlleleCount + 1) - (alleleIndex + refPenalty));
                }
            }


            // approximately rank rankedNewAltAlleleGroup alleles based on sample rankings, if sampleCount ==1 this
            // has already been done
            if (sampleCount > 1)
            {
                std::vector<int> support;
                const unsigned topAllelesCount(rankedNewAltAlleleGroup.size());
                for (unsigned alleleIndex(0); alleleIndex < topAllelesCount; alleleIndex++)
                {
                    const IndelKey& indelKey(rankedNewAltAlleleGroup.key(alleleIndex));
                    const auto indelKeyIter(newAltAlleleKeyScore.find(indelKey));
                    assert(indelKeyIter != newAltAlleleKeyScore.end());
                    support.push_back(indelKeyIter->second);
                }

                const OrthogonalVariantAlleleCandidateGroup tmpCopy(rankedNewAltAlleleGroup);
                rankedNewAltAlleleGroup.clear();
                for (const unsigned altAlleleIndex : sortIndices(support))
                {
                    rankedNewAltAlleleGroup.addVariantAllele(tmpCopy.iter(altAlleleIndex));
                }
            }
        }

        newAltAlleleGroup.clear();

        for (const auto& rankedNewAltAlleleIter : rankedNewAltAlleleGroup.alleles)
        {
            bool isGroupOrthogonal(true);
            for (const auto& newAltAlleleIter : newAltAlleleGroup.alleles)
            {
                if (not is_indel_conflict(rankedNewAltAlleleIter->first, newAltAlleleIter->first))
                {
                    isGroupOrthogonal = false;
                    break;
                }
            }
            if (isGroupOrthogonal)
            {
                newAltAlleleGroup.alleles.push_back(rankedNewAltAlleleIter);
            }
        }
    }

    // put all qualifying alts back together with variants to form an extended allele set:
    OrthogonalVariantAlleleCandidateGroup extendedVariantAlleleGroup(alleleGroup);
    for (const auto& newAltAlleleIter : newAltAlleleGroup.alleles)
    {
        extendedVariantAlleleGroup.alleles.push_back(newAltAlleleIter);
    }

    // rerank and reselect top N alleles, N=callerPloidy per sample, over all samples
    //
    selectTopOrthogonalAllelesInAllSamples(
        sampleCount, callerPloidy, extendedVariantAlleleGroup, alleleGroup, topVariantAlleleIndexPerSample);

    return isEveryAltOrthogonal;
}
