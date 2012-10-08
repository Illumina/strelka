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

#ifndef __STARLING_ALIGN_LIMIT_H
#define __STARLING_ALIGN_LIMIT_H

#include <vector>


// Object which helps to control realignment complexity:
//
struct starling_align_limit {

    starling_align_limit(const unsigned max_alignment_count);

    unsigned
    get_max_toggle(const unsigned n_usable_indels) const {
        if(n_usable_indels>=_max_toggle.size()) {
            return 1;
        } else {
            return _max_toggle[n_usable_indels];
        }
    }

private:
    std::vector<unsigned> _max_toggle;
};


#endif

