// -*- mode: c++; indent-tabs-mode: nil; -*-
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

///
/// \author Sangtae Kim
///

#include <boost/algorithm/string.hpp>
#include <assembly/IterativeAssembler.hh>
#include <blt_util/math_util.hh>
#include "ActiveRegion.hh"

#include "blt_util/fisher_exact_test.hh"

// compile with this macro to get verbose output:
//#define DEBUG_ACTIVE_REGION

void ActiveRegion::processHaplotypes()
{
    // adjust the window size if not enough reads are available
    if (_readBuffer.getEndPos() < _posRange.end_pos)
        _posRange.set_end_pos(_readBuffer.getEndPos());

    for (unsigned sampleId(0); sampleId<_sampleCount; ++sampleId)
    {
        bool isHaplotypingSuccess = processHaplotypesWithCounting(sampleId);
        if (not isHaplotypingSuccess)
        {
            // counting failed. Try assembly.
            isHaplotypingSuccess = processHaplotypesWithAssembly(sampleId);
        }

        if (not isHaplotypingSuccess)
        {
            // both counting and assembly failed
            // do not use haplotyping to determine indel candidacy
            doNotUseHaplotyping();
        }
    }
}

bool ActiveRegion::processHaplotypesWithCounting(unsigned sampleId)
{
    ReadInfo readInfo;
    _readBuffer.getReadSegments(_posRange, readInfo, false);

    unsigned numReads(readInfo.numReads);
    unsigned numReadsCoveringFullRegion((unsigned int) readInfo.readSegments.size());

    // if there are not enough reads fully covering the region, give up counting
    if (numReadsCoveringFullRegion < MinFracReadsCoveringRegion*numReads)
        return false;

    HaplotypeToAlignIdSet haplotypeToAlignIdSet;
    for (const auto& entry : readInfo.readSegments)
    {
        align_id_t alignId = entry.first;
        unsigned currentSampleId = _readBuffer.getSampleId(alignId);
        if (currentSampleId != sampleId) continue;

        const std::string& haplotype(entry.second);

        if (!haplotypeToAlignIdSet.count(haplotype))
            haplotypeToAlignIdSet[haplotype] = std::vector<align_id_t>();
        haplotypeToAlignIdSet[haplotype].push_back(alignId);
    }

#ifdef DEBUG_ACTIVE_REGION
    std::string refStr;
    _ref.get_substring(_posRange.begin_pos, _posRange.size(), refStr);
    std::cerr << _posRange.begin_pos+1 << '\t' << _posRange.end_pos << '\t' << refStr << "\tCounting"<< std::endl;
#endif

    return processSelectedHaplotypes(sampleId, haplotypeToAlignIdSet);
}


bool ActiveRegion::processHaplotypesWithAssembly(unsigned sampleId)
{
    // if reference span is too large, give up assembly
    if (_posRange.size() > MaxRefSpanToBypassAssembly)
    {
        return false;   // assembly fail; bypass indels
    }

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

    std::string refStr;
    _ref.get_substring(_posRange.begin_pos, _posRange.size(), refStr);

    // prefix anchor ends with the first base of the active region
    std::string prefixAnchor;
    _ref.get_substring(beginPos, _posRange.begin_pos - beginPos + 1, prefixAnchor);

    // suffix anchor starts with the last base of the active region
    std::string suffixAnchor;
    _ref.get_substring(_posRange.end_pos-1, endPos-_posRange.end_pos + 1, suffixAnchor);

    // get read segments
    ReadInfo readInfo;
    _readBuffer.getReadSegments(pos_range(beginPos, endPos), readInfo, true);

    // too many reads; do not perform assembly (too time-consuming)
    if (readInfo.numReads > MinNumReadsToBypassAssembly)
        return false;   // assembly fail; bypass indels later

#ifdef DEBUG_ACTIVE_REGION
    std::cerr << _posRange.begin_pos+1 << '\t' << _posRange.end_pos << '\t' << refStr << "\tAssembly"<< std::endl;
#endif
    AssemblyReadInput reads;
    std::vector<align_id_t> readIndexToAlignId;

    for (const auto& entry : readInfo.readSegments)
    {
        align_id_t alignId = entry.first;
        unsigned currentSampleId = _readBuffer.getSampleId(alignId);
        if (currentSampleId != sampleId) continue;

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
    unsigned minWordLength((unsigned int) (prefixAnchor.size() + suffixAnchor.size()));
    assembleOption.minWordLength = minWordLength;

    // maxWordLength must not be smaller than minWordLength.
    unsigned maxWordLength(std::max(minWordLength, ActiveRegion::MaxAssemblyWordSize));
    assembleOption.maxWordLength = maxWordLength;
    assembleOption.minCoverage = MinAssemblyCoverage;

    // perform assembly
    runIterativeAssembler(assembleOption, reads, assemblyReadOutput, contigs);

    HaplotypeToAlignIdSet haplotypeToAlignIdSet;
    unsigned maxHaplotypeLength(0);
    for (unsigned i(0); i<contigs.size(); ++i)
    {
        const std::string& contig(contigs[i].seq);

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

        if (haplotype.length() > maxHaplotypeLength)
            maxHaplotypeLength = (unsigned int) haplotype.length();
        haplotypeToAlignIdSet[haplotype] = std::vector<align_id_t>();
        //unsigned numPseudoReads(0);
        for (unsigned readIndex : contigs[i].supportReads)
        {
            if (assemblyReadOutput[readIndex].isPseudo)
            {
                //++numPseudoReads;
                // TODO: how to add align id for pseudo read?
                // for pseudo reads, assign fake align id
                //haplotypeToAlignIdSet[haplotype].push_back(10000+numPseudoReads);
                continue;
            }

            haplotypeToAlignIdSet[haplotype].push_back(readIndexToAlignId[readIndex]);
        }
    }

    if (haplotypeToAlignIdSet.empty())
        return false;    // assembly fail; bypass indels

    return processSelectedHaplotypes(sampleId, haplotypeToAlignIdSet);
}

void ActiveRegion::doNotUseHaplotyping()
{
#ifdef DEBUG_ACTIVE_REGION
    std::cerr << _posRange.begin_pos+1 << '\t' << _posRange.end_pos << "\tBypass"<< std::endl;
#endif

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



bool ActiveRegion::processSelectedHaplotypes(unsigned sampleId, HaplotypeToAlignIdSet& haplotypeToAlignIdSet)
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
        bool isReference = (haplotype == refStr);
        unsigned count = (unsigned)entry.second.size();

        // count should be no less than MinHaplotypeCount
        if (count < MinHaplotypeCount) continue;
        if (count < secondLargestCount) continue;

        if (count > largestCount)
        {
            secondLargestCount = largestCount;
            largestCount = count;

            secondBestHaplotypePtrList.clear();
            if (bestHaplotypePtr != nullptr)
                secondBestHaplotypePtrList.push_back(std::move(bestHaplotypePtr));

            if (not isReference)
            {
                bestHaplotypePtr = std::unique_ptr<std::string>(new std::string(haplotype));
            }
        }
        else if (count > secondLargestCount)
        {
            secondLargestCount = count;
            secondBestHaplotypePtrList.clear();
            if (not isReference)
            {
                secondBestHaplotypePtrList.push_back(std::unique_ptr<std::string>(new std::string(haplotype)));
            }
        }
        else
        {
            // tie at secondLargestCount
            if (not isReference)
            {
                secondBestHaplotypePtrList.push_back(std::unique_ptr<std::string>(new std::string(haplotype)));
            }
        }
    }

    // now that top two haplotypes are selected we can run various haplotype noise filtration routines:
    //
    if ((largestCount>0) and (secondLargestCount>0))
    {
        // (1) get haplotype 1 and 2 strings:
        const std::string& hap1(bestHaplotypePtr ? *bestHaplotypePtr : refStr);

        std::string* hap2Ptr(nullptr);
        if (not secondBestHaplotypePtrList.empty())
        {
            hap2Ptr = secondBestHaplotypePtrList.front().get();
        }
        const std::string& hap2(hap2Ptr != nullptr ? *hap2Ptr : refStr);

        // (2) identify align-id counts shared in common between the two haplotypes
        const auto hap1MapIter(haplotypeToAlignIdSet.find(hap1));
        assert(hap1MapIter != haplotypeToAlignIdSet.end());
        const auto& hap1AlignIdList(hap1MapIter->second);

        const auto hap2MapIter(haplotypeToAlignIdSet.find(hap2));
        assert(hap2MapIter != haplotypeToAlignIdSet.end());
        const auto& hap2AlignIdList(hap2MapIter->second);

        // build the dup set based on the assumption that the two align_id vectors are sorted:
        std::set<align_id_t> dups;
        {
            auto hap1iter(hap1AlignIdList.begin());
            const auto hap1end(hap1AlignIdList.end());
            auto hap2iter(hap2AlignIdList.begin());
            const auto hap2end(hap2AlignIdList.end());
            while((hap1iter != hap1end) and (hap2iter != hap2end))
            {
                if (*hap1iter < *hap2iter)
                {
                    hap1iter++;
                }
                else if (*hap1iter > *hap2iter)
                {
                    hap2iter++;
                }
                else
                {
                    dups.insert(*hap1iter);
                    hap1iter++;
                    hap2iter++;
                }
            }
        }

        // get new totals counts based support for hapX given that we only consider hap1 and hap2 as possibilities:
        const unsigned hap1UniqueCount(largestCount-dups.size());
        const unsigned hap2UniqueCount(secondLargestCount-dups.size());

#if 0
        const unsigned hap12UniqueCount(hap1UniqueCount+hap2UniqueCount);

        // trial a haplotype evidence ratio filter -- this should be helpful at high depth:
        if (hap12UniqueCount > 10)
        {
            if (safeFrac(hap2UniqueCount,hap12UniqueCount) <= 0.10)
            {
                secondBestHaplotypePtrList.clear();
            }
        }
#endif

        //
        // lambda used to get stranded read counts in downstream noise filtration functions:
        //
        auto getHaplotypeNonDupFwdCount = [&](const std::string& haplotype) {
                unsigned fwdCount(0);
                const auto hapMapIter(haplotypeToAlignIdSet.find(haplotype));
                assert(hapMapIter != haplotypeToAlignIdSet.end());
                const auto& alignIdList(hapMapIter->second);

                for (const auto alignId : alignIdList)
                {
                    if (dups.count(alignId) > 0) continue;
                    if (_readBuffer.getAlignInfo(alignId).isForwardStrand) fwdCount++;
                }
                return fwdCount;
            };

#if 0
        //
        // this is a more general strand bias test assuming the stranded counts are available, it doesn't appear to be specific enough to help
        //
        if (not secondBestHaplotypePtrList.empty())
        {
            // get stranded read counts:
            const unsigned hap1UniqueFwdCount(getHaplotypeNonDupFwdCount(hap1));
            assert(hap1UniqueFwdCount <= hap1UniqueCount);

            const unsigned hap2UniqueFwdCount(getHaplotypeNonDupFwdCount(hap2));
            assert(hap2UniqueFwdCount <= hap2UniqueCount);
            const double biasPval = fisher_exact_test_pval_2x2(hap1UniqueFwdCount, hap2UniqueFwdCount, (hap1UniqueCount-hap1UniqueFwdCount), (hap2UniqueCount-hap2UniqueFwdCount));

            if (biasPval < 0.001)
            {
                // determine which haplotype is the "problem"
                const double hap1AlleleBias=std::abs(0.5-(hap1UniqueFwdCount/ static_cast<double>(hap1UniqueCount)));
                const double hap2AlleleBias=std::abs(0.5-(hap2UniqueFwdCount/ static_cast<double>(hap2UniqueCount)));
                if (hap1AlleleBias > hap2AlleleBias)
                {
                    bestHaplotypePtr.release();
                }
                else
                {
                    secondBestHaplotypePtrList.clear();
                }
            }
        }
#endif

        // test specifically for a very clean sequencer phasing error appearing as haplotype2 when the
        // sample contains a single true homozygous haplotpye which we've detected in haplotype1
        //
        // The error is found by testing for:
        //
        // 1. the top two haplotypes are the same length and differ by only one basecall
        // 2. the second haplotype is observed exclusively on a single strand.
        // 3. that basecall difference is found at the (begin|end) of a homopolymer track, changing the base to match
        //    the hpol base in the hap containing onlyh (rev|fwd) support.
        // 5. the corresponding hpol is at least minPhaseErrorHpolSize in length (this is the hpol in hap1, not reference)
        //
        if (not secondBestHaplotypePtrList.empty())
        {
            static const int minPhaseErrorHpolSize(10);

            // test condition 1:
            if (doHaplotypesMeetPhasingErrorCondition1(hap1,hap2))
            {
                // get stranded read counts:
                const unsigned hap1UniqueFwdCount(getHaplotypeNonDupFwdCount(hap1));
                assert(hap1UniqueFwdCount <= hap1UniqueCount);

                const unsigned hap2UniqueFwdCount(getHaplotypeNonDupFwdCount(hap2));
                assert(hap2UniqueFwdCount <= hap2UniqueCount);

                // test condition 2.:
                if ((hap2UniqueFwdCount == 0) or (hap2UniqueFwdCount == hap2UniqueCount))
                {
                    const auto retval = std::mismatch(hap1.begin(), hap1.end(), hap2.begin());

                    // test conditions 3 and 4:
                    bool isPhaseError(false);
                    if (hap2UniqueFwdCount == 0)
                    {
                        auto hap2iter(retval.second);
                        const char hap2base(*hap2iter);
                        for (; hap2iter != hap2.end(); hap2iter++)
                        {
                            if (*hap2iter != hap2base) break;
                        }
                        if ((hap2iter - retval.second) > minPhaseErrorHpolSize) isPhaseError = true;
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
                        if ((retval.second - hap2iter) > minPhaseErrorHpolSize) isPhaseError = true;

                    }

                    if (isPhaseError)
                    {
                        secondBestHaplotypePtrList.clear();
                    }
                }
            }
        }
    }


    // two haplotypes are selected. Now process them.
    uint8_t haplotypeId(1);
    if (bestHaplotypePtr != nullptr)
    {
        const auto& haplotype(*bestHaplotypePtr);
        const auto& alignIdList(haplotypeToAlignIdSet[haplotype]);
        convertToPrimitiveAlleles(sampleId, haplotype, alignIdList, haplotypeId);
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
            convertToPrimitiveAlleles(sampleId, haplotype, alignIdList, haplotypeId);
            ++haplotypeId;
#ifdef DEBUG_ACTIVE_REGION
            std::cerr << haplotype << '\t' << alignIdList.size() << std::endl;
#endif
        }
    }

    return true;
}

void ActiveRegion::convertToPrimitiveAlleles(
    const unsigned sampleId,
    const std::string& haploptypeSeq,
    const std::vector<align_id_t>& alignIdList,
    const uint8_t haplotypeId)
{
    std::string reference;
    _ref.get_substring(_posRange.begin_pos, _posRange.size(), reference);
    if (reference == haploptypeSeq)
        return;

    pos_t referencePos;
    AlignmentResult<int> result;
    referencePos = _posRange.begin_pos;
    _aligner.align(haploptypeSeq.cbegin(),haploptypeSeq.cend(),reference.cbegin(),reference.cend(),result);

    const ALIGNPATH::path_t& alignPath = result.align.apath;

    pos_t haplotypePosOffset = 0;
    if (result.align.beginPos > 0)
    {
        assert(false && "Unexpected alignment segment");
    }

    unsigned numVariants(0);
    for (unsigned pathIndex(0); pathIndex<alignPath.size(); ++pathIndex)
    {
        const ALIGNPATH::path_segment& pathSegment(alignPath[pathIndex]);
        unsigned segmentLength = pathSegment.length;

        std::unique_ptr<IndelKey> indelKeyPtr = nullptr;
        switch (pathSegment.type)
        {
        case ALIGNPATH::SEQ_MATCH:
            referencePos += segmentLength;
            haplotypePosOffset += segmentLength;
            break;
        case ALIGNPATH::SEQ_MISMATCH:
            for (unsigned i(0); i<segmentLength; ++i)
            {
                if (not _polySites[sampleId].isKeyPresent(referencePos))
                    _polySites[sampleId].getRef(referencePos) = 0;

                addBaseId(haplotypeId, haploptypeSeq[haplotypePosOffset], _polySites[sampleId].getRef(referencePos));

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
                    indelKeyPtr = std::unique_ptr<IndelKey>(new IndelKey(insertPos, INDEL::INDEL, 0, insertSeq.c_str()));
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
                    indelKeyPtr = std::unique_ptr<IndelKey>(new IndelKey(deletePos, INDEL::INDEL, segmentLength));
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

        if (indelKeyPtr != nullptr)
        {
            for (auto alignId : alignIdList)
            {
                IndelObservationData indelObservationData;
                const auto& alignInfo(_readBuffer.getAlignInfo(alignId));
                indelObservationData.iat = alignInfo.indelAlignType;
                indelObservationData.id = alignId;
                _indelBuffer.addIndelObservation(alignInfo.sampleId, {*indelKeyPtr, indelObservationData});
            }
            auto* indelDataPtr(_indelBuffer.getIndelDataPtr(*indelKeyPtr));
            assert(indelDataPtr != nullptr && "Missing indelData");

            // determine whether this indel is candidate or private
            indelDataPtr->isConfirmedInActiveRegion = true;

            indelDataPtr->getSampleData(sampleId).haplotypeId += haplotypeId;

            // TODO: perform candidacy test here
        }
    }
}
