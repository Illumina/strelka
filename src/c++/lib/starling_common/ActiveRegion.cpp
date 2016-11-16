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

void ActiveRegion::insertHaplotypeBase(align_id_t alignId, pos_t pos, const std::string& base)
{
    if (!_alignIdToHaplotype.count(alignId))
    {
        // first occurrence of this alignment
        _alignIdToHaplotype[alignId] = std::string();
        for (int i=_posRange.begin_pos; i<pos; ++i)
            _alignIdToHaplotype[alignId] += missingPrefix;
    }
    _alignIdToHaplotype[alignId] += base;
    if (pos == (_posRange.end_pos-1))
        _alignIdReachingEnd.insert(alignId);
}

void ActiveRegion::processHaplotypes(IndelBuffer& indelBuffer, RangeSet& polySites) const
{
    for (unsigned sampleId(0); sampleId<_sampleCount; ++sampleId)
    {
        // separately process haplotypes per sample
        processHaplotypes(indelBuffer, polySites, sampleId);
    }
}

void ActiveRegion::processHaplotypes(IndelBuffer& indelBuffer, RangeSet& polySites, unsigned sampleId) const
{
    std::map<std::string, std::vector<align_id_t>> haplotypeToAlignIdSet;
    std::vector<std::pair<std::string, align_id_t>> softClippedReads;
    for (const auto& entry : _alignIdToHaplotype)
    {
        align_id_t alignId = entry.first;
        unsigned currentSampleId = _alignIdToAlignInfo[alignId % MaxDepth].sampleId;

        if (currentSampleId != sampleId) continue;

        const std::string& haplotype(entry.second);

        // ignore if the read does not cover the start of the active region
        if (haplotype.empty() || haplotype[0] == missingPrefix) continue;

        // ignore if the read does not reach the end of the active region
        if (_alignIdReachingEnd.find(alignId) == _alignIdReachingEnd.end()) continue;

        // separate soft-clipped reads
        if (_alignIdSoftClipped.find(alignId) != _alignIdSoftClipped.end())
        {
            softClippedReads.push_back(std::pair<std::string, align_id_t>(haplotype, alignId));
            continue;
        }

        if (!haplotypeToAlignIdSet.count(haplotype))
            haplotypeToAlignIdSet[haplotype] = std::vector<align_id_t>();
        haplotypeToAlignIdSet[haplotype].push_back(alignId);
    }

    // match soft-clipped reads to haplotypes
    for (auto& entry : haplotypeToAlignIdSet)
    {
        const std::string& haplotype(entry.first);
        auto& alignIdList(entry.second);

        for (const auto& softClipEntry : softClippedReads)
        {
            const std::string& softClippedRead(softClipEntry.first);
            // checks if the haplotype matches a prefix or suffix
            if (boost::starts_with(softClippedRead, haplotype)
                or boost::ends_with(softClippedRead, haplotype))
            {
                align_id_t alignId(softClipEntry.second);
                alignIdList.push_back(alignId);
            }
        }
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

//    std::cout << "***Sample " << sampleId << std::endl;
//    std::cout << '>' << _posRange.begin_pos+1 << '\t' << _posRange.end_pos << '\t' << _refSeq << '\t' << totalCount << std::endl;
    for (const auto& entry : haplotypeToAlignIdSet)
    {
        const std::string& haplotype(entry.first);
        if (haplotype.empty() || haplotype[0] == missingPrefix) continue;

        // ignore if haplotype is a long homopolymer
        if (haplotype.length() > MaxSNVHpolSize and isHomoPolymer(haplotype)) continue;

        const auto& alignIdList(entry.second);
        auto count = alignIdList.size();

//        if (count >= thirdLargestCount)
//            std::cout << haplotype << '\t' << count << std::endl;
        if (count >= thirdLargestCount and haplotype != _refSeq)
        {
            convertToPrimitiveAlleles(sampleId, haplotype, alignIdList, totalCount, count >= secondLargestCount,
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
    IndelBuffer& indelBuffer,
    RangeSet& polySites) const
{
    AlignmentResult<int> result;
    _aligner.align(haploptypeSeq.begin(),haploptypeSeq.end(),_refSeq.begin(),_refSeq.end(),result);
    const ALIGNPATH::path_t& alignPath = result.align.apath;

    pos_t referencePos = _posRange.begin_pos;
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

                indelKeyPtr = std::unique_ptr<IndelKey>(new IndelKey(insertPos, INDEL::INDEL, 0, insertSeq.c_str()));
                ++numVariants;
                isIndelExist = true;
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

                indelKeyPtr = std::unique_ptr<IndelKey>(new IndelKey(deletePos, INDEL::INDEL, segmentLength));
                ++numVariants;
                isIndelExist = true;
            }
            referencePos += segmentLength;
            break;
        }
        default:
            assert(false && "Unexpected alignment segment");
        }

        if (indelKeyPtr != nullptr)
        {
            auto* indelDataPtr = indelBuffer.getIndelDataPtr(*indelKeyPtr);
            if (indelDataPtr == nullptr)
            {
                // novel indel
                for (auto alignId : alignIdList)
                {
                    IndelObservationData indelObservationData;
                    auto& alignInfo(_alignIdToAlignInfo[alignId % MaxDepth]);
                    indelObservationData.iat = alignInfo.indelAlignType;
                    indelObservationData.id = alignId;
                    indelBuffer.addIndelObservation(alignInfo.sampleId, {*indelKeyPtr, indelObservationData});
                }
                indelDataPtr = indelBuffer.getIndelDataPtr(*indelKeyPtr);
            }
            assert(indelDataPtr != nullptr && "Missing indelData");

            // determine whether this indel is candidate or private
            indelDataPtr->isConfirmedInActiveRegion = true;
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
