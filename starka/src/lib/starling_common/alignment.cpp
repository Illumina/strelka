// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
//

/// \file
///
/// \author Chris Saunders
///

#include "alignment_util.hh"

#include "starling_common/alignment.hh"
#include "starling_common/starling_shared.hh"

#include <iostream>



bool
alignment::
is_overmax(const unsigned max_indel_size) const {

    // check that alignment span and its indels can be handled
    const known_pos_range pr(get_alignment_range(*this));

    // test if total set of alignments in read exceeds
    // read_length+max_indel_size
    //
    if((pr.end_pos-pr.begin_pos)>
       static_cast<pos_t>(apath_read_length(path)+max_indel_size)){
        return true;
    }

    // test if any individual indel exceeds max_indel_size
    using namespace ALIGNPATH;
    const unsigned as(path.size());
    for(unsigned i(0);i<as;++i){
        const path_segment& ps(path[i]);
        if((i==0) || ((i+1)==as)) continue;
        if(((ps.type==INSERT) || (ps.type==DELETE)) &&
           (ps.length>max_indel_size)){
            return true;
        }
    }
    return false;
}



std::ostream&
operator<<(std::ostream& os,
           const alignment& al) {

    os << "ALIGNMENT pos: " << al.pos
       << " strand: " << (al.is_fwd_strand? 'F' : 'R')
       << " path: " << apath_to_cigar(al.path);
    //    if(al.is_overmax(max_indel_size)) os << " overmax";
    if(al.is_seq_swap()) os << " seq_swap";
    os << "\n";

    return os;
}
