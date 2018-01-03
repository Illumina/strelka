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

#include "RegionProcessor.hh"

#include <iostream>



void
RegionProcessor::
addToRegion(
    const std::string& chrom,
    const pos_t outputPos)
{
    if (nullptr == _osptr) return;

    // determine if we need to flush current range:
    if (_is_range)
    {
        if ((chrom != _chrom) || ((_prange.end_pos+1) != outputPos))
        {
            flush();
        }
    }

    if (_is_range)
    {
        _prange.set_end_pos(outputPos);
    }
    else
    {
        _chrom=chrom;
        _prange.set_begin_pos(outputPos-1);
        _prange.set_end_pos(outputPos);
        _is_range=true;
    }
}



void
RegionProcessor::
flush()
{
    if (nullptr == _osptr) return;
    if (! _is_range) return;

    (*_osptr) << _chrom << '\t' << _prange.begin_pos << '\t' << _prange.end_pos << '\n';

    _is_range=false;
}
