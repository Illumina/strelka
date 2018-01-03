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

#include "alignment_util.hh"

#include "starling_common/alignment.hh"

#include <iostream>



bool
alignment::
is_overmax(const unsigned max_indel_size) const
{
    // test if any individual indel exceeds maxIndelSize
    using namespace ALIGNPATH;
    const unsigned as(path.size());
    for (unsigned i(0); i<as; ++i)
    {
        const path_segment& ps(path[i]);
        if ((i==0) || ((i+1)==as)) continue;
        if (((ps.type==INSERT) || (ps.type==DELETE)) &&
            (ps.length>max_indel_size))
        {
            return true;
        }
    }
    return false;
}



std::ostream&
operator<<(std::ostream& os,
           const alignment& al)
{
    os << "ALIGNMENT pos: " << al.pos
       << " strand: " << (al.is_fwd_strand? 'F' : 'R')
       << " path: " << apath_to_cigar(al.path);
    //    if(al.is_overmax(maxIndelSize)) os << " overmax";
    if (al.is_seq_swap()) os << " seq_swap";
    os << "\n";

    return os;
}
