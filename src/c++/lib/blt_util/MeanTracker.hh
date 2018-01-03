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

#pragma once


/// simple online mean utility
struct MeanTracker
{
    double
    mean() const
    {
        return ((_count==0) ? 0. : (_total/_count));
    }

    unsigned
    count() const
    {
        return _count;
    }

    void
    addObservation(const double obs)
    {
        _total += obs;
        _count++;
    }

    void
    clear()
    {
        _total = 0.;
        _count = 0;
    }

private:
    double _total = 0.;
    unsigned _count = 0;
};
