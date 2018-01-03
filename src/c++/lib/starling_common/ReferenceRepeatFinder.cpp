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

#include "ReferenceRepeatFinder.hh"

void ReferenceRepeatFinder::updateRepeatSpan(pos_t pos)
{
    // calculate repeat counter
    auto base(_ref.get_base(pos));
    auto posIndex(pos % _maxBufferSize);

    _isAnchor[posIndex] = true;
    for (auto repeatUnitLength(1u); repeatUnitLength<=_maxRepeatUnitLength; ++repeatUnitLength)
    {
        auto prevBase(_ref.get_base(pos - repeatUnitLength));
        auto repeatUnitIndex(repeatUnitLength-1);
        unsigned repeatSpan;
        if (prevBase != 'N' and base == prevBase)
            repeatSpan = _repeatSpan[(pos - 1) % _maxBufferSize][repeatUnitIndex] + 1;
        else
            repeatSpan = repeatUnitLength;
        _repeatSpan[posIndex][repeatUnitIndex] = repeatSpan;
        if (repeatSpan >= repeatUnitLength*2u and repeatSpan >= _minRepeatSpan)
        {
            // repeat found
            if (repeatSpan == repeatUnitLength*2u or repeatSpan == _minRepeatSpan)
            {
                // unset _isAnchor for pos [pos-repeatSpan+1, pos-1]
                for (pos_t prevPos(pos-1u); prevPos > (pos_t)(pos-repeatSpan); --prevPos)
                {
                    pos_t prevPosIndex(prevPos % _maxBufferSize);
                    _isAnchor[prevPosIndex] = false;
                }
            }
            _isAnchor[posIndex] = false;
        }
    }
}

void ReferenceRepeatFinder::initRepeatSpan(pos_t pos)
{
    pos_t minPos = pos - 2*_maxRepeatUnitLength + 1;
    if (minPos < _ref.get_offset())
        minPos = _ref.get_offset();

    auto posIndex(minPos % _maxBufferSize);

    // initialize _repeatSpan for minPos
    for (auto repeatUnitLength(1u); repeatUnitLength<=_maxRepeatUnitLength; ++repeatUnitLength)
    {
        auto repeatUnitIndex(repeatUnitLength-1);
        _repeatSpan[posIndex][repeatUnitIndex] = repeatUnitLength;
    }

    // Update repeatSpan information up to pos + _maxRepeatUnitLength*2 - 1
    for (pos_t initPos(minPos); initPos < (pos_t)(pos + _maxRepeatUnitLength*2u); ++initPos)
    {
        updateRepeatSpan(initPos);
    }
}
