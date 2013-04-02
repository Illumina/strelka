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

#ifndef __GROUPER_CONTIG_HH
#define __GROUPER_CONTIG_HH

#include "starling_common/alignment.hh"


struct grouper_contig : public alignment {

    grouper_contig() :
        is_usable(true) {}

    void
    clear() {
        alignment::clear();
        id.clear();
        chrom.clear();
        seq.clear();
        is_usable=true;
    }

    std::string id;
    std::string chrom;
    std::string seq;
    bool is_usable;  // < mark contig as usable if it exceeds MAX_CONTIG_SIZE
};


//debugging output:
std::ostream& operator<<(std::ostream& os, const grouper_contig& ctg);

#endif
