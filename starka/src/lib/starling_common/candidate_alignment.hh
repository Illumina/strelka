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

#ifndef __CANDIDATE_ALIGNMENT_HH
#define __CANDIDATE_ALIGNMENT_HH

#include "starling_common/alignment.hh"
#include "starling_common/indel_set.hh"


struct candidate_alignment {

    bool
    operator<(const candidate_alignment& rhs) const {
        if(al < rhs.al) return true;
        if(al == rhs.al) {
            if(leading_indel_key < rhs.leading_indel_key) return true;
            if(leading_indel_key == rhs.leading_indel_key) {
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


#endif
