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

#include "blt_util/blt_types.hh"
#include "blt_util/pos_range.hh"

#include <iosfwd>
#include <string>


/// manage creation of a region track written out in BED format
struct RegionProcessor
{
    explicit
    RegionProcessor(
        std::ostream* osptr) :
        _osptr(osptr)
    {}

    ~RegionProcessor()
    {
        flush();
    }

    /// positions must be added in order:
    void
    addToRegion(
        const std::string& chrom,
        const pos_t outputPos);

    // write out any pending ranges:
    void
    flush();

private:
    std::ostream* _osptr;

    bool _is_range = false;
    std::string _chrom;
    pos_range _prange;
};
