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

#include "OrthogonalVariantAlleleCandidateGroupUtil.hh"
#include "indel_util.hh"
#include "blt_util/prob_util.hh"
#include "blt_util/sort_util.hh"

#include <vector>


#ifdef DEBUG_INDEL_OVERLAP
#include "blt_util/log.hh"

#include <iostream>
#endif


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
    const bool isTier1Only,
    const unsigned minDistanceFromReadEdge)
{
    const unsigned nonrefAlleleCount(alleleGroup.size());
    std::map<unsigned,unsigned> countReadIds;
    for (unsigned nonrefAlleleIndex(0); nonrefAlleleIndex<nonrefAlleleCount; nonrefAlleleIndex++)
    {
        const IndelSampleData& isd(alleleGroup.data(nonrefAlleleIndex).getSampleData(sampleIndex));
        for (const auto& score : isd.read_path_lnp)
        {
            if (isTier1Only && (! score.second.is_tier1_read)) continue;

            if (minDistanceFromReadEdge > 0)
            {
                if (score.second.read_pos < 0) continue;

                // read_pos is the position adjacent to the indel, add one to read_pos from this point to express the
                // concept that the indel is "within" a certain distance from the edge.
                const pos_t intersect_read_pos(score.second.read_pos+1);
                if (intersect_read_pos < static_cast<pos_t>(minDistanceFromReadEdge)) continue;
                if ((score.second.read_length - (intersect_read_pos + 1)) < static_cast<pos_t>(minDistanceFromReadEdge)) continue;
            }

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
    getAlleleGroupUnionReadIds(sampleIndex, alleleGroup, readIds, isTier1Only);
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

    // reference allele is fixed at index 0 by convention
    static const unsigned refAlleleIndex(0);

    // Note that alelleLogLhood values are never initialized - the loops below should be guaranteed
    // to set every alleleLogLhood value, so this is not required.
    alleleLogLhood.resize(fullAlleleCount);

    bool isZeroAlleleCoverage(true);
    bool isPartialAlleleCoverage(false);
    for (unsigned nonrefAlleleIndex(0); nonrefAlleleIndex<nonrefAlleleCount; nonrefAlleleIndex++)
    {
        const auto& indelData(alleleGroup.data(nonrefAlleleIndex));
        const IndelSampleData& indelSampleData(indelData.getSampleData(sampleIndex));

        const auto iditer(indelSampleData.read_path_lnp.find(readId));
        if (iditer==indelSampleData.read_path_lnp.end())
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
        alleleLogLhood[nonrefAlleleIndex+1] = path_lnp.indel;

        isZeroAlleleCoverage=false;
    }

    assert(not isZeroAlleleCoverage);

    // Handle a read which only supports a subset of alleles
    if (isPartialAlleleCoverage)
    {
        for (unsigned nonrefAlleleIndex(0); nonrefAlleleIndex < nonrefAlleleCount; nonrefAlleleIndex++)
        {
            const auto& indelData(alleleGroup.data(nonrefAlleleIndex));
            const IndelSampleData& indelSampleData(indelData.getSampleData(sampleIndex));

            const auto iditer(indelSampleData.read_path_lnp.find(readId));
            if (iditer != indelSampleData.read_path_lnp.end()) continue;

            alleleLogLhood[nonrefAlleleIndex+1] = alleleLogLhood[refAlleleIndex];
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
    normalizeLogDistro(alleleProb.begin(), alleleProb.end(), maxIndex);
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
    static const unsigned refAlleleIndex(0);

    // For each allele, sum the posterior support from all reads
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
        assert(fullAlleleIndex>0);
        const auto& alleleIter(alleleGroup.iter(fullAlleleIndex-1));
        if (not isReferenceRankFound) referenceRank++;
        rankedAlleleGroup.addVariantAllele(alleleIter);
    }

    assert(rankedAlleleGroup.size() == nonrefAlleleCount);
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

    if (topSize<topAlleleGroup.size())
    {
        topAlleleGroup.alleles.resize(topSize);
    }
}



void
selectTopOrthogonalAllelesInAllSamples(
    const unsigned sampleCount,
    const std::vector<unsigned>& callerPloidyPerSample,
    const OrthogonalVariantAlleleCandidateGroup& inputAlleleGroup,
    OrthogonalVariantAlleleCandidateGroup& topAlleleGroup,
    std::vector<unsigned>& topVariantAlleleIndexPerSample)
{
    assert(sampleCount == callerPloidyPerSample.size());
    assert(! inputAlleleGroup.empty());

    topAlleleGroup.clear();

    // This structure is used to create an approximate allele rank over all samples
    std::map<IndelKey,unsigned> topVariantAlleleKeyScore;

    // Record the most likely non-reference allele at this locus in each sample
    std::vector<IndelKey> topVariantAllelePerSample;

    for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
    {
        // In each sample, rank and select top N alleles, N=callerPloidyPerSample, and accumulate these top alleles
        // over all samples in topAlleleGroup
        //
        const unsigned sampleCallerPloidy(callerPloidyPerSample[sampleIndex]);
        OrthogonalVariantAlleleCandidateGroup topAlleleGroupInSample;
        assert(sampleCallerPloidy > 0);
        selectTopOrthogonalAllelesInSample(sampleIndex, inputAlleleGroup, sampleCallerPloidy,
                                           topAlleleGroupInSample);

        // Iterate over the top alleles from this sample to update
        // topAlleleGroup and topVariantAlleleKeyScore
        const unsigned topAllelesInSampleCount(topAlleleGroupInSample.size());
        for (unsigned alleleIndex(0); alleleIndex < topAllelesInSampleCount; alleleIndex++)
        {
            const IndelKey& indelKey(topAlleleGroupInSample.key(alleleIndex));
            auto indelKeyIter(topVariantAlleleKeyScore.find(indelKey));
            if (indelKeyIter == topVariantAlleleKeyScore.end())
            {
                // When allele is observed the first time, set up a new entry in the scoring map, and add the
                // allele to topAlleleGroup
                auto retVal = topVariantAlleleKeyScore.insert(std::make_pair(indelKey,0));
                indelKeyIter = retVal.first;
                topAlleleGroup.addVariantAllele(topAlleleGroupInSample.iter(alleleIndex));
            }

            // Add the (ploidy-adjusted) rank of the allele to the allele's score:
            assert(alleleIndex < sampleCallerPloidy);
            indelKeyIter->second += (sampleCallerPloidy-alleleIndex);
        }

        // Set topVariantAllelePerSample:
        if (topAllelesInSampleCount>0)
        {
            const IndelKey& indelKey(topAlleleGroupInSample.key(0));
            topVariantAllelePerSample.push_back(indelKey);
        }
        else
        {
            topVariantAllelePerSample.push_back(IndelKey::noIndel());
        }
    }

    assert(topVariantAllelePerSample.size() == sampleCount);

    // Rank topAlleleGroup alleles based on the approximate pan-sample support score in topVariantAlleleKeyScore
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

    // Convert top variant per sample from indel key to topAlleleGroup index:
    {
        std::map<IndelKey,unsigned> indelIndex;
        const unsigned topAllelesCount(topAlleleGroup.size());
        for (unsigned alleleIndex(0); alleleIndex < topAllelesCount; alleleIndex++)
        {
            const IndelKey& indelKey(topAlleleGroup.key(alleleIndex));
            indelIndex[indelKey] = alleleIndex;
        }

        topVariantAlleleIndexPerSample.resize(sampleCount);
        std::fill(std::begin(topVariantAlleleIndexPerSample), std::end(topVariantAlleleIndexPerSample), 0);
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            const auto& topIndelKey(topVariantAllelePerSample[sampleIndex]);
            if (topIndelKey == IndelKey::noIndel()) continue;
            const auto indelIter(indelIndex.find(topIndelKey));
            assert(indelIter != indelIndex.end());
            topVariantAlleleIndexPerSample[sampleIndex] = indelIter->second;
        }
    }
}



/// helper for reference_contig_segment
static
void
append_ref_subseq(
    const reference_contig_segment& ref,
    const pos_t start_pos,
    const pos_t end_pos,
    std::string& out_seq)
{
    for (pos_t p(start_pos); p<end_pos; ++p)
    {
        out_seq += ref.get_base(p);
    }
}



/// get the equiv of VCF ALT sequences for each member of an allele group as if these alleles were reported in
/// a single locus record, then mark any repeated alts for filtration
///
/// \return true if any repeats are found
static
bool
getAlleleGroupAltRepeats(
    const reference_contig_segment& ref,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    std::vector<bool>& isRepeatAlt)
{
    bool isFirst(true);
    known_pos_range2 alleleSetRange;
    const unsigned alleleSize(alleleGroup.size());
    for (unsigned alleleIndex(0); alleleIndex < alleleSize; alleleIndex++)
    {
        const IndelKey& indelKey(alleleGroup.key(alleleIndex));
        assert(not indelKey.is_breakpoint());

        if (isFirst)
        {
            alleleSetRange.set_range(indelKey.pos, indelKey.right_pos());
            isFirst = false;
        }
        else
        {
            alleleSetRange.merge_range(known_pos_range2(indelKey.pos, indelKey.right_pos()));
        }
    }

    bool isAnyRepeats(false);
    isRepeatAlt.resize(alleleSize);
    std::set<std::string> alts;
    for (unsigned alleleIndex(0); alleleIndex < alleleSize; alleleIndex++)
    {
        const IndelKey& indelKey(alleleGroup.key(alleleIndex));

        std::string alt;
        append_ref_subseq(ref, alleleSetRange.begin_pos(), indelKey.pos, alt);
        alt += indelKey.insert_seq();
        append_ref_subseq(ref, indelKey.right_pos(), alleleSetRange.end_pos(), alt);

        isRepeatAlt[alleleIndex] = (alts.count(alt)!=0);
        if (not isRepeatAlt[alleleIndex])
        {
            alts.insert(alt);
        }
        else
        {
            isAnyRepeats = true;
        }
    }

    return isAnyRepeats;
}



bool
addAllelesAtOtherPositions(
    const reference_contig_segment& ref,
    const unsigned sampleCount,
    const std::vector<unsigned>& callerPloidy,
    const pos_t pos,
    const pos_t largestTotalIndelRefSpanPerRead,
    const IndelBuffer& indelBuffer,
    OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    std::vector<unsigned>& topVariantAlleleIndexPerSample)
{
    const pos_t minIndelBufferPos(pos-largestTotalIndelRefSpanPerRead);

    const unsigned inputAlleleCount(alleleGroup.size());
    assert(inputAlleleCount!=0);

    bool isEveryAltOrthogonal(true);

    // First get the set of new candidate alt alleles from positions other than 'pos', and store these in
    // newAltAlleleGroup
    //
    OrthogonalVariantAlleleCandidateGroup newAltAlleleGroup;
    {
        const known_pos_range inputAlleleGroupRange(alleleGroup.getReferenceRange());

        // extend end_pos by one to ensure that we find indels adjacent to the right end of the inputAlleleGroupRange
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

            if (altAlleleKey.isMismatch()) continue;

            // must be a candidate allele:
            const IndelData& altAlleleData(getIndelData(altAlleleIter));
            if (not indelBuffer.isCandidateIndel(altAlleleKey, altAlleleData)) continue;

            // Remove indels with do-not-genotype status, these will be written to the output via the forced
            // indel output pathway instead
            if (altAlleleData.doNotGenotype) continue;

            // must be orthogonal to all input alleles:
            /// TODO: comment bug? it looks like this is not a 'must' criteria anymore.
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

#ifdef DEBUG_INDEL_OVERLAP
    log_os << "ZEBRA pos/extended-region-candidate-alleles: " << pos << " " << newAltAlleleGroup << "\n";
#endif

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
                // Rank alt alleles and include from highest to lowest unless interference clique is broken:
                //
                // TODO: This ranking is currently done wrt only the newAltAllele set. It seems like it should be done
                // TODO: wrt newAltAllele+inputAlleleGroup instead.
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
                        // Initialize indel in score map and global allele list if it hasn't been observed
                        // in another sample already:
                        auto retVal = newAltAlleleKeyScore.insert(std::make_pair(indelKey, 0));
                        indelKeyIter = retVal.first;
                        rankedNewAltAlleleGroup.addVariantAllele(rankedNewAltAlleleGroupPerSample.iter(alleleIndex));

                    }
                    if (referenceRank == alleleIndex) refPenalty = 1;
                    indelKeyIter->second += ((newAltAlleleCount + 1) - (alleleIndex + refPenalty));
                }
            }


            // Rank rankedNewAltAlleleGroup alleles based on approximate global scoring computed above in
            // newAltAlleleKeyScore, if sampleCount == 1 the alleles have already been ranked based on read-support.
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

        // Now that the new alleles have been ranked, test them in rank order for conflicting with each other
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

#ifdef DEBUG_INDEL_OVERLAP
        log_os << "ZEBRA pos/ranked-extended-region-candidate-alleles: " << pos << " " << newAltAlleleGroup << "\n";
#endif
    }

    // put all qualifying alts back together with variants to form an extended allele set:
    OrthogonalVariantAlleleCandidateGroup extendedVariantAlleleGroup(alleleGroup);
    for (const auto& newAltAlleleIter : newAltAlleleGroup.alleles)
    {
        extendedVariantAlleleGroup.alleles.push_back(newAltAlleleIter);
    }

#ifdef DEBUG_INDEL_OVERLAP
    log_os << "ZEBRA pos/all-region-candidate-alleles: " << pos << " " << extendedVariantAlleleGroup << "\n";
#endif

    // rerank and reselect top N alleles, N=callerPloidy per sample, over all samples
    //
    selectTopOrthogonalAllelesInAllSamples(
        sampleCount, callerPloidy, extendedVariantAlleleGroup, alleleGroup, topVariantAlleleIndexPerSample);

#ifdef DEBUG_INDEL_OVERLAP
    log_os << "ZEBRA pos/ranked-all-region-candidate-alleles: " << pos << " " << alleleGroup << "\n";
#endif

    const unsigned alleleSize(alleleGroup.size());
    if (alleleSize > 1)
    {
        // check for the rare condition that two alleles would resolve to identical
        // REF/ALT(s) (typically this means that a proximal SNV has not been joined
        // to the candidate indel allele
        //
        // TODO should we print a warning here?
        std::vector<bool> isRepeatAlt;
        if (getAlleleGroupAltRepeats(ref, alleleGroup, isRepeatAlt))
        {
            OrthogonalVariantAlleleCandidateGroup filteredAlleleGroupCopy;
            for (unsigned alleleIndex(0); alleleIndex < alleleSize; alleleIndex++)
            {
                if (isRepeatAlt[alleleIndex]) continue;
                filteredAlleleGroupCopy.alleles.push_back(alleleGroup.iter(alleleIndex));
            }

            // filtration of some alleles requires that we rerank and recompute top-per-sample:
            selectTopOrthogonalAllelesInAllSamples(
                sampleCount, callerPloidy, filteredAlleleGroupCopy, alleleGroup, topVariantAlleleIndexPerSample);
        }
    }

    return isEveryAltOrthogonal;
}
