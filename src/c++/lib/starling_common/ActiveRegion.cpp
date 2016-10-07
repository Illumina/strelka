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
/// \author Sangtae Kim
///

#include <boost/algorithm/string.hpp>
#include <assembly/SmallAssembler.hh>
#include <assembly/IterativeAssembler.hh>
#include "ActiveRegion.hh"

// adoptation of get_snp_hpol_size in blt_common
static unsigned getHomoPolymerSize(const std::string& haplotype, const pos_t pos)
{
    // count upstream repeats:
    bool isUpRepeat(false);
    char upRepeat('N');
    unsigned upSize(0);
    for (pos_t i(pos-1); i>=0; i--)
    {
        if (isUpRepeat)
        {
            if (upRepeat != haplotype[i]) break;
        }
        else
        {
            upRepeat = haplotype[i];
            isUpRepeat = true;
            if (upRepeat == 'N') break;
        }
        upSize++;
    }

    // count downstream repeats:
    bool isDownRepeat(false);
    char downRepeat('N');
    unsigned downSize(0);
    const pos_t haplotypeLength(haplotype.length());
    for (pos_t i(pos+1); i<haplotypeLength; i++)
    {
        if (isDownRepeat)
        {
            if (downRepeat != haplotype[i]) break;
        }
        else
        {
            downRepeat = haplotype[i];
            isDownRepeat = true;
            if (downRepeat == 'N') break;
        }
        downSize++;
    }

    return 1+((downRepeat==upRepeat) ? upSize+downSize : std::max(upSize,downSize) );
}

static bool isHomoPolymer(const std::string& haplotype)
{
    if (haplotype.length() == 0) return true;
    char firstBase = haplotype[0];
    for (unsigned i(1); i<haplotype.length(); ++i)
        if (haplotype[i] != firstBase)
            return false;
    return true;
}

void ActiveRegion::processHaplotypes(IndelBuffer &indelBuffer, RangeSet &polySites) const
{
    for (unsigned sampleId(0); sampleId<_sampleCount; ++sampleId)
    {
        bool isAssemblyRequired = processHaplotypesWithCounting(indelBuffer, polySites, sampleId);
        if (isAssemblyRequired)
            processHaplotypesWithAssembly(indelBuffer, polySites, sampleId);
    }
}

bool ActiveRegion::processHaplotypesWithCounting(IndelBuffer& indelBuffer, RangeSet& polySites, unsigned sampleId) const
{
    HaplotypeInfo haplotypeInfo;
    _readBuffer.getHaplotypeReads(_posRange, haplotypeInfo, false);

    unsigned numReads(haplotypeInfo.numReads);
    unsigned numReadsCoveringFullRegion(haplotypeInfo.readSegments.size());
    unsigned maxHaplotypeLength(0);
    std::map<std::string, std::vector<align_id_t>> haplotypeToAlignIdSet;
    for (const auto& entry : haplotypeInfo.readSegments)
    {
        align_id_t alignId = entry.first;
        unsigned currentSampleId = _readBuffer.getSampleId(alignId);
        if (currentSampleId != sampleId) continue;

        const std::string& haplotype(entry.second);
        if (haplotype.length() > maxHaplotypeLength)
            maxHaplotypeLength = haplotype.length();

        if (!haplotypeToAlignIdSet.count(haplotype))
            haplotypeToAlignIdSet[haplotype] = std::vector<align_id_t>();
        haplotypeToAlignIdSet[haplotype].push_back(alignId);
    }

    bool isAssemblyRequired(false);
    if (std::max(maxHaplotypeLength, _posRange.size()) > 10 and haplotypeInfo.numReads < 1000u and numReadsCoveringFullRegion < 0.5*numReads)
    {
        isAssemblyRequired = true;

        return isAssemblyRequired;
    }

    // determine threshold to select 3 haplotypes with the largest counts
    unsigned largestCount = MinHaplotypeCount;
    unsigned secondLargestCount = MinHaplotypeCount;
    unsigned thirdLargestCount = MinHaplotypeCount;
    unsigned totalCount = 0;
    for (const auto& entry : haplotypeToAlignIdSet)
    {
        auto count = entry.second.size();

        totalCount += count;
        if (count > thirdLargestCount)
        {
            if (count > secondLargestCount)
            {
                if (count > largestCount)
                {
                    thirdLargestCount = secondLargestCount;
                    secondLargestCount = largestCount;
                    largestCount = (unsigned)count;
                }
                else
                {
                    thirdLargestCount = secondLargestCount;
                    secondLargestCount = (unsigned)count;
                }
            }
            else
                thirdLargestCount = (unsigned)count;
        }
    }

    std::string refStr;
    _ref.get_substring(_posRange.begin_pos, _posRange.size(), refStr);
//    std::cout << "chr20" << '\t' << _posRange.begin_pos+1 << '\t' << _posRange.end_pos << '\t' << refStr << "\tCounting"<< std::endl;
    for (const auto& entry : haplotypeToAlignIdSet)
    {
        const std::string& haplotype(entry.first);

        // ignore if haplotype is a long homopolymer
        if (haplotype.length() > MaxSNVHpolSize and haplotype.length() == refStr.length() and isHomoPolymer(haplotype)) continue;

        const auto& alignIdList(entry.second);
        auto count = alignIdList.size();

//        if (count >= secondLargestCount)
//            std::cout << haplotype << '\t' << count << std::endl;

        if (count >= secondLargestCount and haplotype != refStr)
        {
            convertToPrimitiveAlleles(sampleId, haplotype, alignIdList, totalCount, count >= secondLargestCount, _posRange,
                                      indelBuffer, polySites);
        }
    }

    return isAssemblyRequired;
}

void ActiveRegion::processHaplotypesWithAssembly(IndelBuffer& indelBuffer, RangeSet& polySites, unsigned sampleId) const
{
    pos_t beginPos = std::max(_readBuffer.getBeginPos(), (pos_t)(_posRange.begin_pos - AssemblyFlankingSequenceLength));
    pos_t endPos = std::min(_readBuffer.getEndPos(), (pos_t)(_posRange.end_pos + AssemblyFlankingSequenceLength));

    std::string refStr;
    _ref.get_substring(_posRange.begin_pos, _posRange.size(), refStr);

    std::string prefixAnchor;
    _ref.get_substring(beginPos, (_posRange.begin_pos-beginPos), prefixAnchor);

    std::string suffixAnchor;
    _ref.get_substring(_posRange.end_pos, (endPos-_posRange.end_pos), suffixAnchor);

    HaplotypeInfo haplotypeInfo;
    _readBuffer.getHaplotypeReads(pos_range(beginPos, endPos), haplotypeInfo, true);

    if (haplotypeInfo.numReads > 1000u)
        return;

//    std::cout << "chr20" << '\t' << _posRange.begin_pos+1 << '\t' << _posRange.end_pos << '\t' << refStr << "\tAssembly"<< std::endl;
    AssemblyReadInput reads;
    std::vector<align_id_t> readIndexToAlignId;

    for (const auto& entry : haplotypeInfo.readSegments)
    {
        align_id_t alignId = entry.first;
        unsigned currentSampleId = _readBuffer.getSampleId(alignId);
        if (currentSampleId != sampleId) continue;

        const std::string& haplotype(entry.second);
        if (not haplotype.empty())
        {
            reads.push_back(haplotype);
            readIndexToAlignId.push_back(alignId);
        }
    }

    AssemblyReadOutput readInfo;
    Assembly contigs;

    IterativeAssemblerOptions assembleOption;
    assembleOption.minWordLength = 25;
    assembleOption.maxWordLength = 76;
    assembleOption.minCoverage = 3;
    runIterativeAssembler(assembleOption, reads, readInfo, contigs);

    std::map<std::string, std::vector<align_id_t>> haplotypeToAlignIdSet;
    unsigned maxHaplotypeLength(0);
    for (unsigned i(0); i<contigs.size(); ++i)
    {
        const std::string& haplotype(contigs[i].seq);
        if (haplotype.length() > maxHaplotypeLength)
            maxHaplotypeLength = haplotype.length();
        haplotypeToAlignIdSet[haplotype] = std::vector<align_id_t>();
        unsigned numPseudoReads(0);
        for (unsigned readIndex : contigs[i].supportReads)
        {
            if (readInfo[readIndex].isPseudo)
            {
                // TODO: how to add align id for pseudo read?
                ++numPseudoReads;
                haplotypeToAlignIdSet[haplotype].push_back(10000+numPseudoReads);
                continue;
            }

            haplotypeToAlignIdSet[haplotype].push_back(readIndexToAlignId[readIndex]);
        }
    }

    // determine threshold to select 3 haplotypes with the largest counts
    unsigned largestCount = (unsigned)(reads.size()*0.1);
    unsigned secondLargestCount = largestCount;
    unsigned thirdLargestCount = largestCount;
    unsigned totalCount = 0;
    for (const auto& entry : haplotypeToAlignIdSet)
    {
        auto count(entry.second.size());
        totalCount += count;
        if (count > thirdLargestCount)
        {
            if (count > secondLargestCount)
            {
                if (count > largestCount)
                {
                    thirdLargestCount = secondLargestCount;
                    secondLargestCount = largestCount;
                    largestCount = (unsigned)count;
                }
                else
                {
                    thirdLargestCount = secondLargestCount;
                    secondLargestCount = (unsigned)count;
                }
            }
            else
                thirdLargestCount = (unsigned)count;
        }
    }

    for (const auto& entry : haplotypeToAlignIdSet)
    {
        const std::string& contig(entry.first);

        auto start(contig.find(prefixAnchor));
        if (start == std::string::npos) continue;

        start += prefixAnchor.length();
        auto end(contig.rfind(suffixAnchor));
        if (end == std::string::npos or start >= end) continue;

        const std::string haplotype(contig.substr(start, end-start));

        // ignore if haplotype is a long homopolymer
        if (haplotype.length() > MaxSNVHpolSize and isHomoPolymer(haplotype)) continue;

        const auto& alignIdList(entry.second);
        auto count(alignIdList.size());
//        if (count >= secondLargestCount)
//            std::cout << haplotype << '\t' << count << std::endl;

        if (count >= secondLargestCount)
        {
            convertToPrimitiveAlleles(sampleId, haplotype, alignIdList, totalCount, count >= secondLargestCount, _posRange,
                                      indelBuffer, polySites);
        }
    }
}


void ActiveRegion::convertToPrimitiveAlleles(
    const unsigned sampleId,
    const std::string& haploptypeSeq,
    const std::vector<align_id_t>& alignIdList,
    const unsigned totalReadCount,
    const bool isTopTwo,
    const pos_range posRange,
    IndelBuffer& indelBuffer,
    RangeSet& polySites) const
{
    std::string reference;
    _ref.get_substring(posRange.begin_pos, posRange.size(), reference);
    if (reference == haploptypeSeq)
        return;

    pos_t referencePos;
    AlignmentResult<int> result;
    referencePos = posRange.begin_pos;
    _aligner.align(haploptypeSeq.cbegin(),haploptypeSeq.cend(),reference.cbegin(),reference.cend(),result);

    const ALIGNPATH::path_t& alignPath = result.align.apath;

    pos_t haplotypePosOffset = 0;
    if (result.align.beginPos > 0)
    {
        assert(false && "Unexpected alignment segment");
    }

    std::vector<pos_t> mismatchPositions;
    std::vector<pos_t> mismatchHaplotypePositions;
    unsigned numVariants(0);
    bool isIndelExist(false);
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
                if (_posRange.is_pos_intersect(referencePos))
                    mismatchPositions.push_back(referencePos);
                mismatchHaplotypePositions.push_back(haplotypePosOffset);
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
                auto insertSeq(haploptypeSeq.substr(haplotypePosOffset, segmentLength));
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
                    isIndelExist = true;
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
                    isIndelExist = true;
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
            auto* indelDataPtr = indelBuffer.getIndelDataPtr(*indelKeyPtr);
//            if (indelDataPtr == nullptr)
            {
                // novel indel
                for (auto alignId : alignIdList)
                {
                    IndelObservationData indelObservationData;
                    const auto& alignInfo(_readBuffer.getAlignInfo(alignId));
                    indelObservationData.iat = alignInfo.indelAlignType;
                    indelObservationData.id = alignId;
                    indelBuffer.addIndelObservation(alignInfo.sampleId, {*indelKeyPtr, indelObservationData});
                }
                indelDataPtr = indelBuffer.getIndelDataPtr(*indelKeyPtr);
            }
            assert(indelDataPtr != nullptr && "Missing indelData");

            // determine whether this indel is candidate or private
            indelDataPtr->isConfirmedInActiveRegion = true;
//            std::cout << indelKeyPtr->pos << "," << indelKeyPtr->insertSequence
//                      << "," << indelKeyPtr->deletionLength
//                      << "," << indelBuffer.isCandidateIndel(*indelKeyPtr, *indelDataPtr) << std::endl;
            // TODO: perform candidacy test here
        }
    }

    if (!isTopTwo or mismatchPositions.empty()) return;

    // relax MMDF
    auto readCount = alignIdList.size();
    bool isLowDepth = (readCount < HighDepth) and (readCount < (unsigned)(totalReadCount*HaplotypeFrequencyThreshold));
    if (numVariants == 1 and isLowDepth)
        return;

    for (unsigned i(0); i<mismatchPositions.size(); ++i)
    {
        bool isLongHpol = getHomoPolymerSize(haploptypeSeq, mismatchHaplotypePositions[i]) > MaxSNVHpolSize;

        // for SNV only haplotypes, ignore if applying the SNV makes a long homo polymer
        if (!isIndelExist and isLongHpol) continue;

        // register MMDF relax site
        polySites[sampleId].getRef(mismatchPositions[i]) = 1;
    }
}
