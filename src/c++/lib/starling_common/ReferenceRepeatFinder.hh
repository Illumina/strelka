//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
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

/// \file
/// \author Sangtae Kim
///

#pragma once

#include <blt_util/blt_types.hh>
#include <blt_util/reference_contig_segment.hh>
#include <vector>

/// Helper object for ActiveRegionReadBuffer which identifies STR tracks in the reference
///
/// This object's primary role is to find "anchor" positions in the reference, where an anchor
/// means a position of sufficiently high local sequence complexity that it can be used as an
/// edge to an active region.
///
class ReferenceRepeatFinder
{
public:
    ReferenceRepeatFinder(
        const reference_contig_segment& ref,
        const unsigned maxRepeatUnitLength,
        const unsigned maxBufferSize,
        const unsigned minRepeatSpan):
        _ref(ref), _maxRepeatUnitLength(maxRepeatUnitLength),
        _maxBufferSize(maxBufferSize), _minRepeatSpan(minRepeatSpan),
        _repeatSpan(maxBufferSize, std::vector<unsigned>(maxRepeatUnitLength)),
        _isAnchor(maxBufferSize)
    {
    }

    /// Checks if pos is an anchor position (not belonging to a repeat)
    /// \param pos reference position
    /// \return true if pos is an anchor position, false otherwise
    bool isAnchor(pos_t pos) const
    {
        return _isAnchor[pos % _maxBufferSize];
    }

    /// Read upstream and downstream and determine whether pos is within a repeat.
    /// by updating _repeatSpan for [pos-2*_maxRepeatSpan+1,pos+2*_maxRepeatSpan-1]
    void initRepeatSpan(pos_t pos);

    /// Update _repeatSpan for pos assuming _repeatSpan has records for [pos-4*_maxRepeatSpan+2,pos-1]
    /// If a repeat is detected at pos, set _isAnchor false appropriately for previous positions
    void updateRepeatSpan(pos_t pos);

private:
    const reference_contig_segment& _ref;
    const unsigned _maxRepeatUnitLength;
    const unsigned _maxBufferSize;
    const unsigned _minRepeatSpan;

    /// _repeatSpan[posIndex][repeatUnitIndex] records
    /// the longest span of repeats ending at pos (posIndex = pos % _maxBufferSize)
    /// with the repeat unit length repeatUnitIndex (repeatUnitIntex = repeatUnitLength -1)
    /// E.g. for _ref = "TATATACCCCCAATGAAAAA", _maxBufferSize = 20
    /// _repeatSpan[5][1] = 6 (TATATA) and _repeatSpan[5][2]=3 (ATA)
    /// if the _repeatSpan value is no less than repeatUnitLength*2,
    /// the position is within a repeat
    std::vector<std::vector<unsigned>> _repeatSpan;
    std::vector<bool> _isAnchor;
};
