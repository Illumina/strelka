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
#include <assembly/IterativeAssembler.hh>
#include "ActiveRegion.hh"

// compile with this macro to get verbose output:
//#define DEBUG_ACTIVE_REGION

static bool isHomoPolymer(const std::string& haplotype)
{
    if (haplotype.length() == 0) return true;
    char firstBase = haplotype[0];
    for (unsigned i(1); i<haplotype.length(); ++i)
        if (haplotype[i] != firstBase)
            return false;
    return true;
}

void ActiveRegion::processHaplotypes(IndelBuffer& indelBuffer, RangeSet& polySites)
{
    // adjust the window size if not enough reads are available
    if (_readBuffer.getEndPos() < _posRange.end_pos)
        _posRange.set_end_pos(_readBuffer.getEndPos());

    for (unsigned sampleId(0); sampleId<_sampleCount; ++sampleId)
    {
        bool isHaplotypingSuccess = processHaplotypesWithCounting(indelBuffer, polySites, sampleId);
        if (not isHaplotypingSuccess)
        {
            // counting failed. Try assembly.
            isHaplotypingSuccess = processHaplotypesWithAssembly(indelBuffer, polySites, sampleId);
        }

        if (not isHaplotypingSuccess)
        {
            // both counting and assembly failed
            // bypass indels parsed from BAM
            bypassIndelsInBam(indelBuffer, sampleId);
        }
    }
}

bool ActiveRegion::processHaplotypesWithCounting(IndelBuffer& indelBuffer, RangeSet& polySites, unsigned sampleId) const
{
    ReadInfo readInfo;
    _readBuffer.getReadSegments(_posRange, readInfo, false);

    unsigned numReads(readInfo.numReads);
    unsigned numReadsCoveringFullRegion((unsigned int) readInfo.readSegments.size());

    // if the fraction of reads fully covering the region is low, give up counting
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

    // determine threshold to select 2 haplotypes with the largest counts
    unsigned largestCount = MinHaplotypeCount;
    unsigned secondLargestCount = MinHaplotypeCount;
    for (const auto& entry : haplotypeToAlignIdSet)
    {
        auto count = entry.second.size();

        if (count > secondLargestCount)
        {
            if (count > largestCount)
            {
                secondLargestCount = largestCount;
                largestCount = (unsigned)count;
            }
            else
            {
                secondLargestCount = (unsigned)count;
            }
        }
    }

    std::string refStr;
    _ref.get_substring(_posRange.begin_pos, _posRange.size(), refStr);

#ifdef DEBUG_ACTIVE_REGION
    std::cerr << _posRange.begin_pos+1 << '\t' << _posRange.end_pos << '\t' << refStr << "\tCounting"<< std::endl;
#endif

    unsigned numNonRefBestHaps(0);
    unsigned numNonRefSecondBestHaps(0);
    for (const auto& entry : haplotypeToAlignIdSet)
    {
        const std::string& haplotype(entry.first);
        if (haplotype == refStr)
            continue;

        auto count = entry.second.size();
        if (count >= largestCount)
            ++numNonRefBestHaps;
        if (count >= secondLargestCount)
            ++numNonRefSecondBestHaps;
    }

    unsigned countThreshold;
    if (numNonRefSecondBestHaps <= 2)
        countThreshold = secondLargestCount;
    else if (numNonRefBestHaps <= 2)
        countThreshold = largestCount;
    else return true;

    uint8_t haplotypeId(0);
    for (const auto& entry : haplotypeToAlignIdSet)
    {
        const std::string& haplotype(entry.first);

        // ignore if haplotype is a long homopolymer
        if (haplotype.length() > MaxSNVHpolSize and haplotype.length() == refStr.length() and isHomoPolymer(haplotype)) continue;

        const auto& alignIdList(entry.second);
        auto count = alignIdList.size();

#ifdef DEBUG_ACTIVE_REGION
        if (count >= countThreshold)
            std::cerr << haplotype << '\t' << count << std::endl;
#endif

        if (count >= countThreshold and haplotype != refStr)
        {
            convertToPrimitiveAlleles(sampleId, haplotype, alignIdList, haplotypeId++,
                                      indelBuffer, polySites);
        }
    }

    // counting success
    return true;
}

bool ActiveRegion::processHaplotypesWithAssembly(IndelBuffer& indelBuffer, RangeSet& polySites, unsigned sampleId) const
{
    // if reference span is too large, give up assembly
    if (_posRange.size() > MaxRefSpanToPerformAssembly)
    {
#ifdef DEBUG_ACTIVE_REGION
        std::cerr << _posRange.begin_pos+1 << '\t' << _posRange.end_pos << "\tAssembly"<< std::endl;
#endif
        return false;   // assembly fail; bypass indels later
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
        unsigned numPseudoReads(0);
        for (unsigned readIndex : contigs[i].supportReads)
        {
            if (assemblyReadOutput[readIndex].isPseudo)
            {
                ++numPseudoReads;
                // TODO: how to add align id for pseudo read?
                // for pseudo reads, assign fake align id
                haplotypeToAlignIdSet[haplotype].push_back(10000+numPseudoReads);
                continue;
            }

            haplotypeToAlignIdSet[haplotype].push_back(readIndexToAlignId[readIndex]);
        }
    }

    if (haplotypeToAlignIdSet.empty())
        return true;    // assembly fail; do not bypass indels

    // determine threshold to select 2 haplotypes with the largest counts
    unsigned largestCount = MinHaplotypeCount;
    unsigned secondLargestCount = largestCount;
    for (const auto& entry : haplotypeToAlignIdSet)
    {
        auto count(entry.second.size());
        if (count > secondLargestCount)
        {
            if (count > largestCount)
            {
                secondLargestCount = largestCount;
                largestCount = (unsigned)count;
            }
            else
            {
                secondLargestCount = (unsigned)count;
            }
        }
    }

    unsigned numNonRefBestHaps(0);
    unsigned numNonRefSecondBestHaps(0);
    for (const auto& entry : haplotypeToAlignIdSet)
    {
        const std::string& haplotype(entry.first);
        if (haplotype == refStr)
            continue;

        auto count = entry.second.size();
        if (count >= largestCount)
            ++numNonRefBestHaps;
        if (count >= secondLargestCount)
            ++numNonRefSecondBestHaps;
    }

    unsigned countThreshold;
    if (numNonRefSecondBestHaps <= 2)
        countThreshold = secondLargestCount;
    else if (numNonRefBestHaps <= 2)
        countThreshold = largestCount;
    else return true;

    uint8_t haplotypeId(0);
    for (const auto& entry : haplotypeToAlignIdSet)
    {
        const std::string& haplotype(entry.first);

        const auto& alignIdList(entry.second);
        auto count(alignIdList.size());

#ifdef DEBUG_ACTIVE_REGION
        if (count >= countThreshold)
            std::cerr << haplotype << '\t' << count << std::endl;
#endif

        if (count >= countThreshold and haplotype != refStr)
        {
            convertToPrimitiveAlleles(sampleId, haplotype, alignIdList, haplotypeId++,
                                      indelBuffer, polySites);
        }
    }

    return true;    // assembly success
}

// simply bypass all indels in BAM
// The original indel candidacy method will be used.
void ActiveRegion::bypassIndelsInBam(IndelBuffer& indelBuffer, unsigned /*sampleId*/) const
{
    auto it(indelBuffer.positionIterator(_posRange.begin_pos));
    const auto it_end(indelBuffer.positionIterator(_posRange.end_pos));

    for (; it!=it_end; ++it)
    {
        const IndelKey& indelKey(it->first);
        IndelData& indelData(getIndelData(it));

        if (indelKey.is_breakpoint()) continue;
        indelData.isConfirmedInActiveRegion = true;
    }
}

void ActiveRegion::convertToPrimitiveAlleles(
    const unsigned sampleId,
    const std::string& haploptypeSeq,
    const std::vector<align_id_t>& alignIdList,
    const uint8_t haplotypeId,
    IndelBuffer& indelBuffer,
    RangeSet& polySites) const
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
                if (not polySites[sampleId].isKeyPresent(referencePos))
                    polySites[sampleId].getRef(referencePos) = 0;

                HaplotypeId complexAlleleId(haplotypeId+1);
                addBaseId(complexAlleleId, haploptypeSeq[haplotypePosOffset], polySites[sampleId].getRef(referencePos));

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
                indelBuffer.addIndelObservation(alignInfo.sampleId, {*indelKeyPtr, indelObservationData});
            }
            auto* indelDataPtr(indelBuffer.getIndelDataPtr(*indelKeyPtr));
            assert(indelDataPtr != nullptr && "Missing indelData");

            // determine whether this indel is candidate or private
            indelDataPtr->isConfirmedInActiveRegion = true;

            indelDataPtr->getSampleData(sampleId).haplotypeId += (haplotypeId+1);

            // TODO: perform candidacy test here
        }
    }
}
