//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
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
/// \author Sangtae Kim
///

#include "ActiveRegion.hh"

#include "assembly/IterativeAssembler.hh"
#include "blt_util/algo_util.hh"

#include "boost/algorithm/string.hpp"
#include "boost/make_unique.hpp"

// compile with this macro to get verbose output:
//#define DEBUG_ACTIVE_REGION

void ActiveRegion::processHaplotypes()
{
    // Check whether the active region is included in the read buffer
    const bool isRangeValid = (_posRange.begin_pos >= _readBuffer.getBeginPos())
                              && (_posRange.end_pos <= _readBuffer.getEndPos());

    // if the active region is not included in the read buffer or if it is too large,
    // bypass haplotyping
    if ((! isRangeValid) || (_posRange.size() > MaxRefSpanToBypassAssembly))
    {
        doNotUseHaplotyping();
    }
    else
    {
        bool isHaplotypingSuccess = processHaplotypesWithCounting();
        if (not isHaplotypingSuccess)
        {
            // counting failed. Try assembly.
            isHaplotypingSuccess = processHaplotypesWithAssembly();
        }

        if (not isHaplotypingSuccess)
        {
            // both counting and assembly failed
            // do not use haplotyping to determine indel candidacy
            doNotUseHaplotyping();
        }
    }
}

bool ActiveRegion::processHaplotypesWithCounting()
{
    ReadInfo readInfo;
    _readBuffer.getReadSegments(_posRange, readInfo, false);

    unsigned numReads(readInfo.numReads);
    unsigned numReadsCoveringFullRegion((unsigned int) readInfo.readSegments.size());

    // if there are not enough reads fully covering the region, give up counting
    if ((numReads == 0) or (numReadsCoveringFullRegion < MinFracReadsCoveringRegion*numReads))
        return false;

    HaplotypeToAlignIdSet haplotypeToAlignIdSet;
    for (const auto& entry : readInfo.readSegments)
    {
        align_id_t alignId = entry.first;

        const std::string& haplotype(entry.second);

        if (!haplotypeToAlignIdSet.count(haplotype))
            haplotypeToAlignIdSet[haplotype] = std::vector<align_id_t>();
        haplotypeToAlignIdSet[haplotype].push_back(alignId);
    }

#ifdef DEBUG_ACTIVE_REGION
    std::string refStr;
    _ref.get_substring(_posRange.begin_pos, _posRange.size(), refStr);
    std::cerr << _sampleIndex << "\t" << _posRange.begin_pos+1 << '\t' << _posRange.end_pos << '\t' << refStr << "\tCounting"<< std::endl;
#endif

    return processSelectedHaplotypes(haplotypeToAlignIdSet, numReads);
}


bool ActiveRegion::processHaplotypesWithAssembly()
{
    // Expand the region to include left/right anchors.
    // TODO: anchors may be too short if there are SNVs close to anchors
    // prefix anchor
    pos_t minBeginPos(0u);
    if (_posRange.begin_pos > ActiveRegionReadBuffer::MaxAssemblyPadding)
        minBeginPos = _posRange.begin_pos - ActiveRegionReadBuffer::MaxAssemblyPadding;
    minBeginPos = std::max(_readBuffer.getBeginPos(), minBeginPos);

    pos_t beginPos(_posRange.begin_pos);
    for (; beginPos > minBeginPos; --beginPos)
    {
        // anchor should not include a variant position
        if (_readBuffer.isCandidateVariant(beginPos-1)) break;
    }

    // suffix anchor
    pos_t maxEndPos = std::min(_readBuffer.getEndPos(), _posRange.end_pos + ActiveRegionReadBuffer::MaxAssemblyPadding);
    pos_t endPos(_posRange.end_pos);
    for (; endPos < maxEndPos; ++endPos)
    {
        if (_readBuffer.isCandidateVariant(endPos)) break;
    }

    // prefix anchor ends with the first base of the active region
    std::string prefixAnchor;
    _ref.get_substring(beginPos, _posRange.begin_pos - beginPos + 1, prefixAnchor);

    // suffix anchor starts with the last base of the active region
    std::string suffixAnchor;
    _ref.get_substring(_posRange.end_pos-1, endPos-_posRange.end_pos + 1, suffixAnchor);

    unsigned minReadSegmentLength((unsigned int) (prefixAnchor.size() + suffixAnchor.size()));

    // get read segments
    ReadInfo readInfo;
    _readBuffer.getReadSegments(pos_range(beginPos, endPos), readInfo, true, minReadSegmentLength);

    // too many reads; do not perform assembly (too time-consuming)
    if (readInfo.numReads > MinNumReadsToBypassAssembly)
        return false;   // assembly fail; bypass indels later

    AssemblyReadInput reads;
    std::vector<align_id_t> readIndexToAlignId;

    for (const auto& entry : readInfo.readSegments)
    {
        align_id_t alignId = entry.first;

        const std::string& readSegment(entry.second);

        if (not readSegment.empty())
        {
            reads.push_back(readSegment);
            readIndexToAlignId.push_back(alignId);
        }
    }

    AssemblyReadOutput assemblyReadOutput;
    Assembly contigs;

    IterativeAssemblerOptions assembleOption;

    // We only accept haplotypes that start with prefixAnchor and end with suffixAnchor.
    // So, the minimum size of legitimate haplotypes is prefixAnchor.size() + suffixAnchor.size().
    // For most regions, minWordLength will be 20.
    // If there are SNVs closer to anchors, minWordLength may be smaller than 20.
    assembleOption.minWordLength = minReadSegmentLength;

    // maxWordLength must not be smaller than minWordLength.
    unsigned maxWordLength(std::max(minReadSegmentLength, ActiveRegion::MaxAssemblyWordSize));
    assembleOption.maxWordLength = maxWordLength;
    assembleOption.minCoverage = MinAssemblyCoverage;

    // perform assembly
    runIterativeAssembler(assembleOption, reads, assemblyReadOutput, contigs);

    unsigned totalNumReads(0);
    for (auto assemblyReadInfo : assemblyReadOutput)
    {
        if (assemblyReadInfo.isUsed and (not assemblyReadInfo.isPseudo))
            ++totalNumReads;
    }

    HaplotypeToAlignIdSet haplotypeToAlignIdSet;
    bool isNonRefHaplotypeFound(false);

    std::string refStr;
    _ref.get_substring(_posRange.begin_pos, _posRange.size(), refStr);

#ifdef DEBUG_ACTIVE_REGION
    std::cerr << _sampleIndex << "\t" << _posRange.begin_pos+1 << '\t' << _posRange.end_pos << '\t' << refStr << "\tAssembly"<< std::endl;
#endif

    for (unsigned contigIndex(0); contigIndex<contigs.size(); ++contigIndex)
    {
        const std::string& contig(contigs[contigIndex].seq);
        // ignore if the contig does not contain prefix anchor
        auto start(contig.find(prefixAnchor));
        if (start == std::string::npos) continue;

        // remove prefix padding
        start += prefixAnchor.length() - 1;

        // ignore if the contig does not contain suffix anchor
        auto end(contig.rfind(suffixAnchor));
        if (end == std::string::npos or start > end) continue;

        // remove suffix padding
        end += 1;

        const std::string haplotype(contig.substr(start, end-start));

        auto alignIds = std::vector<align_id_t>();
        bool containsUniqueRead(false);
        for (unsigned readIndex : contigs[contigIndex].supportReads)
        {
            const auto& assemblyReadInfo(assemblyReadOutput[readIndex]);
            if (assemblyReadInfo.isPseudo)
            {
                // pseudo reads are ignored
                continue;
            }
            if ((not containsUniqueRead) and (assemblyReadInfo.contigIds.size() == 1))
                containsUniqueRead = true;
            alignIds.push_back(readIndexToAlignId[readIndex]);
        }

        // ignore if there's no read uniquely supporting the contig
        if (not containsUniqueRead) continue;

        if (haplotype != refStr)
            isNonRefHaplotypeFound = true;
        haplotypeToAlignIdSet[haplotype] = alignIds;
    }

    // assembly fails if no alt haplotype is found
    if (not isNonRefHaplotypeFound)
        return false;

    return processSelectedHaplotypes(haplotypeToAlignIdSet, totalNumReads);
}

void ActiveRegion::doNotUseHaplotyping()
{
#ifdef DEBUG_ACTIVE_REGION
    std::cerr << _sampleIndex << "\t" << _posRange.begin_pos+1 << '\t' << _posRange.end_pos << "\tBypass"<< std::endl;
#endif

    assert(_posRange.end_pos > _posRange.begin_pos);

    auto it(_indelBuffer.positionIterator(_posRange.begin_pos));
    const auto it_end(_indelBuffer.positionIterator(_posRange.end_pos));

    for (; it!=it_end; ++it)
    {
        const IndelKey& indelKey(it->first);
        IndelData& indelData(getIndelData(it));

        if (indelKey.is_breakpoint()) continue;
        indelData.isConfirmedInActiveRegion = true;
    }
}



/// \return true if haplotypes are the same length and have exactly one mismatch
static
bool
doHaplotypesMeetPhasingErrorCondition1(
    const std::string& hap1,
    const std::string& hap2)
{
    // 2.
    if ((hap1.length() == hap2.length()) and (hap1 != hap2))
    {
        const auto retval = std::mismatch(hap1.begin(), hap1.end(), hap2.begin());
        const auto retval2 = std::mismatch(retval.first + 1, hap1.end(), retval.second + 1);
        if (retval2.first == hap1.end())
        {
            return true;
        }
    }
    return false;
}


/// \return true if haplotype2 appears to be a sequencer phasing error 'echo' of haplotype 1
///
/// test specifically for a very clean sequencer phasing error appearing as haplotype2 when the
/// sample contains a single true homozygous haplotype which we've detected in haplotype1
///
/// The error is found by testing for:
///
/// 1. the top two haplotypes are the same length and differ by only one basecall
/// 2. the second haplotype is observed exclusively on a single strand.
/// 3. that basecall difference is found at the (begin|end) of a homopolymer track, changing the base to match
///    the hpol base in the hap containing only (rev|fwd) support.
/// 4. the corresponding hpol is at least minPhaseErrorHpolSize in length (this is the hpol in hap1, not reference)
///
static
bool
isFilterSecondHaplotypeAsSequencerPhasingNoise(
    const ActiveRegionReadBuffer& readBuffer,
    const HaplotypeToAlignIdSet& haplotypeToAlignIdSet,
    const std::string& hap1,
    const std::string& hap2)
{
    static const int minPhaseErrorHpolSize(10);

    // test condition 1
    //
    // arrange this filter to come first so that we don't have to identify dups or read strands for most cases:
    if (not doHaplotypesMeetPhasingErrorCondition1(hap1,hap2)) return false;


    // test condition 2
    //

    // 2a: get haplotype supporting read lists:
    const auto hap1MapIter(haplotypeToAlignIdSet.find(hap1));
    assert(hap1MapIter != haplotypeToAlignIdSet.end());
    const auto& hap1AlignIdList(hap1MapIter->second);

    const auto hap2MapIter(haplotypeToAlignIdSet.find(hap2));
    assert(hap2MapIter != haplotypeToAlignIdSet.end());
    const auto& hap2AlignIdList(hap2MapIter->second);

    // 2b: identify duplicate reads:
    const std::set<align_id_t> dups(
        getDuplicatesInSortedInput(std::begin(hap1AlignIdList),std::end(hap1AlignIdList),
                                   std::begin(hap2AlignIdList),std::end(hap2AlignIdList)));

    // 2c: get total and unique counts for haplotype2:
    const unsigned hap2Count(hap2AlignIdList.size());
    const unsigned hap2UniqueCount(hap2Count-dups.size());

    // 2d: identify stranded counts:

    /// \return unique fwd-strand counts
    auto getHaplotypeNonDupFwdCount = [&](const std::vector<align_id_t>& alignIdList)
    {
        unsigned fwdCount(0);
        for (const auto alignId : alignIdList)
        {
            if (dups.count(alignId) > 0) continue;
            if (readBuffer.getAlignInfo(alignId).isForwardStrand) fwdCount++;
        }
        return fwdCount;
    };

    const unsigned hap2UniqueFwdCount(getHaplotypeNonDupFwdCount(hap2AlignIdList));
    assert(hap2UniqueFwdCount <= hap2UniqueCount);

    // 2e: test condition 2
    if ((hap2UniqueFwdCount > 0) and (hap2UniqueFwdCount < hap2UniqueCount)) return false;


    // test conditions 3 and 4
    //

    const auto retval = std::mismatch(hap1.begin(), hap1.end(), hap2.begin());
    if (hap2UniqueFwdCount == 0)
    {
        auto hap2iter(retval.second);
        const char hap2base(*hap2iter);
        for (; hap2iter != hap2.end(); hap2iter++)
        {
            if (*hap2iter != hap2base) break;
        }
        return ((hap2iter - retval.second) > minPhaseErrorHpolSize);
    }
    else
    {
        auto hap2iter(retval.second);
        const char hap2base(*hap2iter);
        while (true)
        {
            if (*hap2iter != hap2base) break;
            if (hap2iter == hap2.begin()) break;
            hap2iter--;
        }
        return ((retval.second - hap2iter) > minPhaseErrorHpolSize);
    }
}



bool ActiveRegion::processSelectedHaplotypes(
    HaplotypeToAlignIdSet& haplotypeToAlignIdSet,
    const unsigned totalNumReads)
{
    // determine threshold to select 2 haplotypes with the largest counts
    unsigned largestCount(0);
    unsigned secondLargestCount(0);

    // haplotype with the most count
    std::unique_ptr<std::string> bestHaplotypePtr;

    // haplotype with the second most count (size >1 if there are ties)
    std::vector<std::unique_ptr<std::string>> secondBestHaplotypePtrList;

    std::string refStr;
    _ref.get_substring(_posRange.begin_pos, _posRange.size(), refStr);

    for (const auto& entry : haplotypeToAlignIdSet)
    {
        const std::string& haplotype(entry.first);
        const bool isReference = (haplotype == refStr);
        unsigned count = (unsigned)entry.second.size();

        // count should be no less than MinHaplotypeCount
        if (count < MinHaplotypeCount) continue;
        if (count < secondLargestCount) continue;

        if (count > largestCount)
        {
            secondLargestCount = largestCount;
            largestCount = count;

            secondBestHaplotypePtrList.clear();
            if (bestHaplotypePtr)
                secondBestHaplotypePtrList.push_back(std::move(bestHaplotypePtr));

            if (not isReference)
            {
                bestHaplotypePtr.reset(new std::string(haplotype));
            }
        }
        else if (count > secondLargestCount)
        {
            secondLargestCount = count;
            secondBestHaplotypePtrList.clear();
            if (not isReference)
            {
                secondBestHaplotypePtrList.push_back(boost::make_unique<std::string>(haplotype));
            }
        }
        else
        {
            // tie at secondLargestCount
            if (not isReference)
            {
                secondBestHaplotypePtrList.push_back(boost::make_unique<std::string>(haplotype));
            }
        }
    }

    // now that top haplotypes are selected we can run haplotype noise filtration routines:
    //
    if ((largestCount>0) and (secondLargestCount>0) and (not secondBestHaplotypePtrList.empty()))
    {
        //
        // get haplotype 1 and 2 strings:
        //
        const std::string& hap1(bestHaplotypePtr ? *bestHaplotypePtr : refStr);

        std::string* hap2Ptr(nullptr);
        if (not secondBestHaplotypePtrList.empty())
        {
            hap2Ptr = secondBestHaplotypePtrList.front().get();
        }
        const std::string& hap2(hap2Ptr != nullptr ? *hap2Ptr : refStr);

        if (isFilterSecondHaplotypeAsSequencerPhasingNoise(_readBuffer, haplotypeToAlignIdSet, hap1, hap2))
        {
            secondBestHaplotypePtrList.clear();
        }
    }

    // top haplotypes are selected. Now process them.
    //
    uint8_t haplotypeId(1);
    if (bestHaplotypePtr)
    {
        const auto& haplotype(*bestHaplotypePtr);
        const auto& alignIdList(haplotypeToAlignIdSet[haplotype]);
        convertToPrimitiveAlleles(haplotype, alignIdList, totalNumReads, haplotypeId);
        ++haplotypeId;
#ifdef DEBUG_ACTIVE_REGION
        std::cerr << haplotype << '\t' << alignIdList.size() << std::endl;
#endif
    }

    if ((not secondBestHaplotypePtrList.empty()) and ((haplotypeId+secondBestHaplotypePtrList.size()-1) <= 2))
    {
        // including second best haplotypes doesn't exceed 2 alt haplotypes
        for (const auto& haplotypePtr : secondBestHaplotypePtrList)
        {
            const auto& haplotype(*haplotypePtr);
            const auto& alignIdList(haplotypeToAlignIdSet[haplotype]);
            convertToPrimitiveAlleles(haplotype, alignIdList, totalNumReads, haplotypeId);
            ++haplotypeId;
#ifdef DEBUG_ACTIVE_REGION
            std::cerr << haplotype << '\t' << alignIdList.size() << std::endl;
#endif
        }
    }

    return true;
}

void ActiveRegion::convertToPrimitiveAlleles(
    const std::string& haploptypeSeq,
    const std::vector<align_id_t>& alignIdList,
    const unsigned totalNumReads,
    const uint8_t haplotypeId)
{
    std::string reference;
    _ref.get_substring(_posRange.begin_pos, _posRange.size(), reference);
    if (reference == haploptypeSeq)
        return;

    pos_t referencePos;
    AlignmentResult<int> result;
    referencePos = _posRange.begin_pos;

    // \TODO: this aligner is already left-shifting, why was the extra left-shift logic added below?
    _aligner.align(haploptypeSeq.cbegin(),haploptypeSeq.cend(),reference.cbegin(),reference.cend(),result);

    const ALIGNPATH::path_t& alignPath = result.align.apath;

    pos_t haplotypePosOffset = 0;
    if (result.align.beginPos > 0)
    {
        assert(false && "Unexpected alignment segment");
    }

    unsigned numVariants(0);
    const float altHaplotypeCountRatio(alignIdList.size()/static_cast<float>(totalNumReads));
    for (const auto& pathSegment : alignPath)
    {
        const unsigned segmentLength = pathSegment.length;

        std::unique_ptr<IndelKey> indelKeyPtr;
        switch (pathSegment.type)
        {
        case ALIGNPATH::SEQ_MATCH:
            referencePos += segmentLength;
            haplotypePosOffset += segmentLength;
            break;
        case ALIGNPATH::SEQ_MISMATCH:
            for (unsigned i(0); i<segmentLength; ++i)
            {
                _candidateSnvBuffer.addCandidateSnv(_sampleIndex, referencePos, haploptypeSeq[haplotypePosOffset], haplotypeId, altHaplotypeCountRatio);

                ++referencePos;
                ++haplotypePosOffset;
            }
            ++numVariants;
            break;
        case ALIGNPATH::INSERT:
        {
            if (segmentLength <= _maxIndelSize)
            {
                // left-align insertion
                // the insertion can be moved left by 1 base
                // if the last base of insertSeq equals to prevBase
                // E.g. GT -> GT(ATAT) vs GT -> G(TATA)T
                pos_t insertPos(referencePos);
                auto insertSeq(haploptypeSeq.substr((unsigned long) haplotypePosOffset, segmentLength));
                char prevBase = _ref.get_base(insertPos-1);
                while (insertSeq.back() == prevBase)
                {
                    // move insertion 1 base to left
                    insertSeq = prevBase + insertSeq;
                    insertSeq.pop_back();

                    --insertPos;
                    prevBase = _ref.get_base(insertPos-1);
                }

                if (prevBase != 'N')
                {
                    indelKeyPtr.reset(new IndelKey(insertPos, INDEL::INDEL, 0, insertSeq.c_str()));
                    ++numVariants;
                }
            }
            haplotypePosOffset += segmentLength;
            break;
        }
        case ALIGNPATH::DELETE:
        {
            if (segmentLength <= _maxIndelSize)
            {
                // left-align deletion
                // the deletion can be moved left by 1 base
                // if the last base of the deleted sequence (lastDeletionBase) equals to prevBase
                // E.g. GT(ATAT) -> GT vs G(TATA)T -> GT
                pos_t deletePos(referencePos);
                char prevBase = _ref.get_base(deletePos-1);
                char lastDeletionBase = _ref.get_base(deletePos + segmentLength - 1);
                while (lastDeletionBase == prevBase)
                {
                    // move deletion 1 base to left
                    --deletePos;
                    lastDeletionBase = _ref.get_base(deletePos + segmentLength - 1);
                    prevBase = _ref.get_base(deletePos-1);
                }

                if (prevBase != 'N')
                {
                    indelKeyPtr.reset(new IndelKey(deletePos, INDEL::INDEL, segmentLength));
                    ++numVariants;
                }
            }
            referencePos += segmentLength;
            break;
        }
        case ALIGNPATH::SOFT_CLIP:
        {
            referencePos += segmentLength;
            break;
        }
        default:
            assert(false && "Unexpected alignment segment");
        }

        if (indelKeyPtr)
        {
            for (const auto alignId : alignIdList)
            {
                IndelObservationData indelObservationData;
                const auto& alignInfo(_readBuffer.getAlignInfo(alignId));
                indelObservationData.iat = alignInfo.indelAlignType;
                indelObservationData.id = alignId;
                _indelBuffer.addIndelObservation(alignInfo.sampleIndex, {*indelKeyPtr, indelObservationData});
            }
            auto* indelDataPtr(_indelBuffer.getIndelDataPtr(*indelKeyPtr));
            assert(indelDataPtr != nullptr && "Missing indelData");

            // determine whether this indel is candidate or private
            indelDataPtr->isConfirmedInActiveRegion = true;

            indelDataPtr->getSampleData(_sampleIndex).haplotypeId += haplotypeId;

            // TODO: why is this a plus? we can't add ratios....
            indelDataPtr->getSampleData(_sampleIndex).altAlleleHaplotypeCountRatio += altHaplotypeCountRatio;

            // TODO: perform candidacy test here
        }
    }
}
