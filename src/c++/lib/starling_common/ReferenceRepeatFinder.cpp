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

#include "ReferenceRepeatFinder.hh"

void ReferenceRepeatFinder::updateRepeatSpan(pos_t pos)
{
    // calculate repeat counter
    auto base(_ref.get_base(pos));
    auto posIndex(pos % MaxBufferSize);

    _isAnchor[posIndex] = true;
    for (auto repeatUnitLength(1u); repeatUnitLength<=MaxRepeatUnitLength; ++repeatUnitLength)
    {
        auto prevBase(_ref.get_base(pos - repeatUnitLength));
        auto repeatUnitIndex(repeatUnitLength-1);
        unsigned repeatSpan;
        if (prevBase != 'N' and base == prevBase)
            repeatSpan = _repeatSpan[(pos - 1) % MaxBufferSize][repeatUnitIndex] + 1;
        else
            repeatSpan = repeatUnitLength;
        _repeatSpan[posIndex][repeatUnitIndex] = repeatSpan;
        if (repeatSpan >= repeatUnitLength*2 and repeatSpan >= MinRepeatSpan)
        {
            // repeat found
            if (repeatSpan == repeatUnitLength*2 or repeatSpan == MinRepeatSpan)
            {
                // unset _isAnchor for pos [pos-repeatCount+1, pos-1]
                for (pos_t prevPos(pos-1u); prevPos > (pos_t)(pos-repeatSpan); --prevPos)
                {
                    pos_t prevPosIndex(prevPos % MaxBufferSize);
                    _isAnchor[prevPosIndex] = false;
                }
            }
            _isAnchor[posIndex] = false;
        }
    }
}

void ReferenceRepeatFinder::initRepeatSpan(pos_t pos)
{
    auto posIndex(pos % MaxBufferSize);

    for (auto repeatUnitLength(1u); repeatUnitLength<=MaxRepeatUnitLength; ++repeatUnitLength)
    {
        auto repeatUnitIndex(repeatUnitLength-1);
        _repeatSpan[posIndex][repeatUnitIndex] = repeatUnitLength;
    }

    for (pos_t initPos(pos); initPos < (pos_t)(pos + MaxRepeatUnitLength*2u); ++initPos)
    {
        updateRepeatSpan(initPos);
    }
}