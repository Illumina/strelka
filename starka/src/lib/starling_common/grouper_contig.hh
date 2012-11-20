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
