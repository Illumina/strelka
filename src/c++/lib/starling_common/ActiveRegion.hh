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

#pragma once

#include "blt_util/blt_types.hh"
#include "starling_common/starling_types.hh"
#include "alignment/GlobalNoClippingAligner.hh"
#include "blt_util/align_path.hh"
#include "IndelBuffer.hh"

#include <string>
#include <map>
#include <set>

/// Represent all haplotypes found in the current active region
class ActiveRegion
{
public:
    ActiveRegion(pos_t start, pos_t end, unsigned numVariants, const std::string& refSeq):
        _start(start), _end(end), _numVariants(numVariants), _refSeq(refSeq),
        _alignIdToHaplotype(),
        _alignIdReachingEnd(),
        _scores(1,-4,-5,-1,-1000),
        _aligner(_scores)
    {}

    pos_t getStart() const
    {
        return _start;
    }
    pos_t getEnd() const
    {
        return _end;
    }
    unsigned getLength() const
    {
        return _end - _start + 1;
    }
    unsigned getNumVariants() const
    {
        return _numVariants;
    }
    bool contains(pos_t pos) const
    {
        return pos >= _start && pos <= _end;
    }
    void insertHaplotypeBase(align_id_t align_id, pos_t pos, const std::string& base);
    void addComplexAllelesToIndelBuffer(IndelBuffer& indelBuffer, std::set<pos_t>& polySites) const;
    const std::string& getReferenceSeq() const { return _refSeq; }

private:
    pos_t _start;
    pos_t _end;
    unsigned _numVariants;
    const std::string& _refSeq;

    std::map<align_id_t, std::string> _alignIdToHaplotype;
    std::set<align_id_t> _alignIdReachingEnd;
    AlignmentScores<int> _scores;
    GlobalNoClippingAligner<int> _aligner;

    void addPrimitiveAllelesToIndelBuffer(
            const std::string &haploptypeSeq,
            std::vector<align_id_t>& alignIdList,
            IndelBuffer& indelBuffer,
            std::set<pos_t>& polySites) const;
};
