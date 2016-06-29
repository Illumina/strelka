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
#include "IndelKey.hh"
#include "indel.hh"

#include <iostream>


void ActiveRegion::insertHaplotypeBase(align_id_t align_id, pos_t pos, const std::string& base)
{
    if (!_alignIdToHaplotype.count(align_id))
    {
        _alignIdToHaplotype[align_id] = std::string();
        for (int i=_start; i<pos; ++i)
            _alignIdToHaplotype[align_id] += '.';
    }
    _alignIdToHaplotype[align_id] += base;
    if (pos == _end)
        _alignIdReachingEnd.insert(align_id);
}

void ActiveRegion::createComplexAlleles() const
{
    std::map<std::string, unsigned> haplotypeCounter;
    unsigned maxCount = 0;
    for (const auto& entry : _alignIdToHaplotype)
    {
        align_id_t alignId = entry.first;
        std::string haplotype = entry.second;
        if (_alignIdReachingEnd.find(alignId) == _alignIdReachingEnd.end())
            haplotype += "*";

        if (!haplotypeCounter.count(haplotype))
            haplotypeCounter[haplotype] = 1;
        else
        {
            unsigned count = haplotypeCounter[haplotype] + 1;
            if (count > maxCount)
                maxCount = count;
            haplotypeCounter[haplotype] = count;
        }
    }

    for (const auto& entry : haplotypeCounter)
    {
        std::string haplotype = entry.first;
        unsigned count = entry.second;
//        if ((count >= 3) && (count >= maxCount/4))
        if (count >= 2)
        {
            if (haplotype.length() == 0 || haplotype[0] == '.' || haplotype[haplotype.length()-1] == '*') continue;
            printVariants(haplotype, _refSeq);
        }
    }
}

void ActiveRegion::printVariants(const std::string &haploptypeSeq, const std::string &referenceSeq) const
{
    AlignmentResult<int> result;
    _aligner.align(haploptypeSeq.begin(),haploptypeSeq.end(),referenceSeq.begin(),referenceSeq.end(),result);
    ALIGNPATH::path_t alignPath = result.align.apath;

//    std::cout << referenceSeq << '\t' << haploptypeSeq << '\t' << alignPath << '\t' << result.align.beginPos << '\t' << result.score << std::endl;
    pos_t pos = _start+1;
    pos_t referenceIndex = 0;
    pos_t haplotypeIndex = 0;
    if (result.align.beginPos > 0)
    {
        std::cout << "DELETE" << '\t' << pos+1 << '\t' << referenceSeq.substr(0, result.align.beginPos) << '\t' << std::endl;
        pos += result.align.beginPos;
        referenceIndex += result.align.beginPos;
    }

    for (unsigned pathIndex(0); pathIndex<alignPath.size(); ++pathIndex)
    {
        const ALIGNPATH::path_segment &pathSegment(alignPath[pathIndex]);
        unsigned segmentLength = pathSegment.length;

        switch (pathSegment.type)
        {
            case ALIGNPATH::SEQ_MATCH:
                pos += segmentLength;
                referenceIndex += segmentLength;
                haplotypeIndex += segmentLength;
                break;
            case ALIGNPATH::SEQ_MISMATCH:
                for (unsigned i(0); i<segmentLength; ++i)
                {
                    std::cout << "MISMATCH" << '\t' << pos++ << '\t' << referenceSeq[referenceIndex++];
                    std::cout << '\t' << haploptypeSeq[haplotypeIndex++] << std::endl;
                }
                break;
            case ALIGNPATH::INSERT:
            case ALIGNPATH::SOFT_CLIP:
                std::cout << "INSERT" << '\t' << pos << '\t' << haploptypeSeq.substr(haplotypeIndex, segmentLength) << std::endl;
                haplotypeIndex += segmentLength;
                break;
            case ALIGNPATH::DELETE:
                std::cout << "DELETE" << '\t' << pos << '\t' << referenceSeq.substr(referenceIndex, segmentLength) << std::endl;
                pos += segmentLength;
                referenceIndex += segmentLength;
                break;
            default:
                std::cout << "Wrong2: " << pathSegment.type << std::endl;
        }
    }
}



