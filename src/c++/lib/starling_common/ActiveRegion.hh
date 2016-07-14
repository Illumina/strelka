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

typedef RangeMap<pos_t,unsigned char> RangeSet;

/// Represent all haplotypes found in the current active region
class ActiveRegion
{
public:
    const float HaplotypeFrequencyThreshold = 0.4; // minimum haplotype frequency to be considered in MMDF relaxation
    const char missingPrefix = '.';

    ActiveRegion(pos_t start, pos_t end, const std::string& refSeq, const GlobalNoClippingAligner<int>& aligner):
        _start(start), _end(end), _refSeq(refSeq),
        _aligner(aligner),
        _alignIdToHaplotype(),
        _alignIdReachingEnd()
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
    bool contains(pos_t pos) const
    {
        return pos >= _start && pos <= _end;
    }
    void insertHaplotypeBase(align_id_t align_id, pos_t pos, const std::string& base);
    void processHaplotypes(IndelBuffer &indelBuffer, RangeSet &polySites) const;
    const std::string& getReferenceSeq() const { return _refSeq; }

private:
    pos_t _start;
    pos_t _end;
    const std::string& _refSeq;
    const GlobalNoClippingAligner<int> _aligner;

    std::map<align_id_t, std::string> _alignIdToHaplotype;
    std::set<align_id_t> _alignIdReachingEnd;

    void convertToPrimitiveAlleles(
            const std::string &haploptypeSeq,
            const std::vector<align_id_t> &alignIdList,
            IndelBuffer &indelBuffer,
            RangeSet &polySites) const;
};
