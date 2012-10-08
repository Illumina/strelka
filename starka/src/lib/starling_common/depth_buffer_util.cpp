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


#include "depth_buffer_util.hh"



void
add_alignment_to_depth_buffer(const alignment& al,
                              depth_buffer& db){

    using namespace ALIGNPATH;

    pos_t ref_head_pos(al.pos);

    const unsigned as(al.path.size());
    for(unsigned i(0);i<as;++i){
        const path_segment& ps(al.path[i]);
        if(ps.type == MATCH) {
            for(unsigned j(0);j<ps.length;++j) db.inc(ref_head_pos+static_cast<pos_t>(j));
        }

        if(is_segment_type_ref_length(ps.type)) ref_head_pos += ps.length;
    }

}

