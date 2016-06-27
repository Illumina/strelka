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
#include "alignment/GlobalAligner.hh"
#include "blt_util/align_path.hh"

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

void ActiveRegion::printHaplotypeSequences() const
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
        if ((count >= 3) && (count >= maxCount/4))
        {
            if (haplotype[0] == '.' || haplotype[haplotype.length()-1] == '*') continue;
            std::cout << haplotype.c_str() << '\t' << count;
            testAlign(haplotype, _refSeq);
            std::cout << std::endl;
        }
        testAlign("AGAGCTCCGGTAGC", "GCTCCGGCAGC");
        std::cout << std::endl;
    }
}

void ActiveRegion::testAlign(const std::string& seq, const std::string& ref) const
{
    AlignmentScores<int> scores(1,-4,-6,-1,-4);
    GlobalAligner<int> aligner(scores);
    AlignmentResult<int> result;
    aligner.align(seq.begin(),seq.end(),ref.begin(),ref.end(),result);
    std::cout << '\t' << apath_to_cigar(result.align.apath);
    std::cout << '\t' << result.align.apath;
}



