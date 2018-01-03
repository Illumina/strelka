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

///
/// \author Chris Saunders
///

#pragma once

#include <cassert>

#include <iosfwd>


struct SamplePloidyState
{
    explicit
    SamplePloidyState(const int ploidy = 2)
    {
        setPloidy(ploidy);
    }

    void
    setPloidy(const int ploidy)
    {
        assert(ploidy<=2 && ploidy>=-1);
        _ploidy = ploidy;
    }

    int
    getPloidy() const
    {
        return _ploidy;
    }

    bool
    isDiploid() const
    {
        return (2 == _ploidy);
    }

    bool
    isHaploid() const
    {
        return (1 == _ploidy);
    }

    bool
    isNoploid() const
    {
        return (0 == _ploidy);
    }

    bool
    isContinuous() const
    {
        return (-1 == _ploidy);
    }

    void
    reset()
    {
        _ploidy = 2;
    }

private:
    /// ploidy: current accept up to 2. 1 is haploid, 0 is "noploid" (inside a deletion) and -1 is continuous
    int _ploidy = 2;
};

std::ostream&
operator<<(std::ostream& os, const SamplePloidyState& sps);
