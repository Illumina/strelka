// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

/// \file
///
/// \author Chris Saunders
///

#include "alignment_util.hh"

#include "starling_common/alignment.hh"
#include "starling_common/starling_base_shared.hh"

#include <iostream>



bool
alignment::
is_overmax(const unsigned max_indel_size) const
{

    // test if any individual indel exceeds max_indel_size
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
    //    if(al.is_overmax(max_indel_size)) os << " overmax";
    if (al.is_seq_swap()) os << " seq_swap";
    os << "\n";

    return os;
}
