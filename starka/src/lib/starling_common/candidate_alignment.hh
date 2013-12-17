// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
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

#pragma once

#include "starling_common/alignment.hh"
#include "starling_common/indel_set.hh"


struct candidate_alignment {

    bool
    operator<(const candidate_alignment& rhs) const {
        if (al < rhs.al) return true;
        if (al == rhs.al) {
            if (leading_indel_key < rhs.leading_indel_key) return true;
            if (leading_indel_key == rhs.leading_indel_key) {
                return (trailing_indel_key < rhs.trailing_indel_key);
            }
        }
        return false;
    }

    alignment al;
    indel_key leading_indel_key;
    indel_key trailing_indel_key;
};


std::ostream& operator<<(std::ostream& os, const candidate_alignment& cal);



// get the keys of the indels present in the candidate alignment
//
void
get_alignment_indels(const candidate_alignment& cal,
                     const unsigned max_indel_size,
                     indel_set_t& indels);
