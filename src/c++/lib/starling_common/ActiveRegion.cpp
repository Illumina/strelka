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

#include "ActiveRegion.hh"

void ActiveRegion::insertHaplotypeBase(align_id_t align_id, pos_t pos, const std::string& base)
{
    if (!_alignIdToHaplotype.count(align_id))
    {
        // first occurrence of this alignment
        _alignIdToHaplotype[align_id] = std::string();
        for (int i=_start; i<pos; ++i)
            _alignIdToHaplotype[align_id] += missingPrefix;
    }
    _alignIdToHaplotype[align_id] += base;
    if (pos == _end)
        _alignIdReachingEnd.insert(align_id);
}

// decompose haplotypes into primitive alleles
void ActiveRegion::processHaplotypes(IndelBuffer& indelBuffer, RangeSet& polySites) const
{
    std::map<std::string, std::vector<align_id_t>> haplotypeToAlignIdSet;
    for (const auto& entry : _alignIdToHaplotype)
    {
        align_id_t alignId = entry.first;
        const std::string& haplotype(entry.second);

        // ignore if the read does not reach the end of the active region
        if (_alignIdReachingEnd.find(alignId) == _alignIdReachingEnd.end()) continue;

        if (!haplotypeToAlignIdSet.count(haplotype))
            haplotypeToAlignIdSet[haplotype] = std::vector<align_id_t>();
        haplotypeToAlignIdSet[haplotype].push_back(alignId);
    }

    // determine threshold to select 2 haplotypes with the largest counts
    unsigned largestCount = 0;
    unsigned secondLargestCount = 0;
    unsigned totalCount = 0;
    for (const auto& entry : haplotypeToAlignIdSet)
    {
        const std::string& haplotype(entry.first);
        if (haplotype.empty() || haplotype[0] == missingPrefix) continue;

        auto count = entry.second.size();
        totalCount += count;
        if (count > secondLargestCount)
        {
            if (count > largestCount)
            {
                secondLargestCount = largestCount;
                largestCount = (unsigned)count;
            }
            else
                secondLargestCount = (unsigned)count;
        }
    }

    for (const auto& entry : haplotypeToAlignIdSet)
    {
        const std::string& haplotype(entry.first);
        if (haplotype.empty() || haplotype[0] == missingPrefix) continue;

        const auto& alignIdList(entry.second);
        auto numReads = alignIdList.size();
        if (numReads >= secondLargestCount && numReads >= totalCount*HaplotypeFrequencyThreshold)
        {
//            std::cout << haplotype << std::endl;
            convertToPrimitiveAlleles(haplotype, alignIdList, indelBuffer, polySites);
        }
    }
}

void ActiveRegion::convertToPrimitiveAlleles(
    const std::string& haploptypeSeq,
    const std::vector<align_id_t>& /*alignIdList*/,
    IndelBuffer& /*indelBuffer*/,
    RangeSet& polySites) const
{
    AlignmentResult<int> result;
    _aligner.align(haploptypeSeq.begin(),haploptypeSeq.end(),_refSeq.begin(),_refSeq.end(),result);
    const ALIGNPATH::path_t& alignPath = result.align.apath;

    pos_t referencePos = _start;
    pos_t haplotypePosOffset = 0;
    if (result.align.beginPos > 0)
    {
//        IndelObservation indelObservation;
//        indelObservation.key.pos = pos;
//        indelObservation.key.type = INDEL::DELETE;
//        indelObservation.key.length = result.align.beginPos;
//        indelObservation.data.is_discovered_in_active_region = true;
//        indelBuffer.addIndelObservation(sampleId, indelObservation);

        referencePos += result.align.beginPos;
    }

    for (unsigned pathIndex(0); pathIndex<alignPath.size(); ++pathIndex)
    {
        const ALIGNPATH::path_segment& pathSegment(alignPath[pathIndex]);
        unsigned segmentLength = pathSegment.length;

        switch (pathSegment.type)
        {
        case ALIGNPATH::SEQ_MATCH:
            referencePos += segmentLength;
            haplotypePosOffset += segmentLength;
            break;
        case ALIGNPATH::SEQ_MISMATCH:
            for (unsigned i(0); i<segmentLength; ++i)
            {
//                    std::cout << "Poly\t" << (pos+1) << std::endl;
                polySites.getRef(referencePos) = 1;
                ++referencePos;
                ++haplotypePosOffset;
            }
            break;
        case ALIGNPATH::INSERT:
        case ALIGNPATH::SOFT_CLIP:
        {
//                IndelObservation indelObservation;
//                indelObservation.key.pos = pos;
//                indelObservation.key.type = INDEL::INSERT;
//                indelObservation.key.length = segmentLength;
//                indelObservation.data.insert_seq = haploptypeSeq.substr(haplotypeIndex, segmentLength);
//                indelObservation.data.is_discovered_in_active_region = true;
//                indelBuffer.addIndelObservation(sampleId, indelObservation);
            haplotypePosOffset += segmentLength;
            break;
        }
        case ALIGNPATH::DELETE:
        {
//                IndelObservation indelObservation;
//                indelObservation.key.pos = pos;
//                indelObservation.key.type = INDEL::DELETE;
//                indelObservation.key.length = segmentLength;
//                indelObservation.data.is_discovered_in_active_region = true;
//                indelBuffer.addIndelObservation(sampleId, indelObservation);

            referencePos += segmentLength;
            break;
        }
        default:
            assert(false && "Unexpected alignment segment");
        }
    }
}



