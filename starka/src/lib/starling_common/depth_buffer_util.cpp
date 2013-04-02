// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

/// \file
///
/// \author Chris Saunders
///


#include "depth_buffer_util.hh"



void
add_alignment_to_depth_buffer(const alignment& al,
                              depth_buffer& db) {

    using namespace ALIGNPATH;

    pos_t ref_head_pos(al.pos);

    const unsigned as(al.path.size());
    for(unsigned i(0); i<as; ++i) {
        const path_segment& ps(al.path[i]);
        if(ps.type == MATCH) {
            for(unsigned j(0); j<ps.length; ++j) db.inc(ref_head_pos+static_cast<pos_t>(j));
        }

        if(is_segment_type_ref_length(ps.type)) ref_head_pos += ps.length;
    }

}

