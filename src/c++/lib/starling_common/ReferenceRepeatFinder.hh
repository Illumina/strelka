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

#include <blt_util/blt_types.hh>
#include <blt_util/reference_contig_segment.hh>

class ReferenceRepeatFinder
{
public:
    static const pos_t MinRepeatSpan = (pos_t)3u;

    // maximum repeat unit to consider
    static const unsigned MaxRepeatUnitLength = 50u;

    // maximum buffer size in bases (must be larger than the maximum read size + max indel size
    static const unsigned MaxBufferSize = 1000u;

    ReferenceRepeatFinder(const reference_contig_segment& ref): _ref(ref) {}

    /// checks if pos is an anchor position
    /// \param pos reference position
    /// \return true if pos is an anchor position, false otherwise
    bool isAnchor(pos_t pos) const
    {
        return _isAnchor[pos % MaxBufferSize];
    }

    void initRepeatSpan(pos_t pos);
    void updateRepeatSpan(pos_t pos);

private:
    const reference_contig_segment& _ref;

    // to record reference repeat count
    unsigned _repeatSpan[MaxBufferSize][MaxRepeatUnitLength];
    bool _isAnchor[MaxBufferSize];
};