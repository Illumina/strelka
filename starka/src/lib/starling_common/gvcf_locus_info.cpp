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


#include "starling_common/gvcf_locus_info.hh"

#include <iostream>


void
shared_modifiers::
write_filters(std::ostream& os) const {

    if(_filters.none()) {
        os << "PASS";
        return;
    }

    bool is_sep(false);
    for(unsigned i(0);i<VCF_FILTERS::SIZE;++i){
        if(_filters[i]) {
            if(is_sep) { os << ";"; }
            else       { is_sep=true; }
            os << VCF_FILTERS::get_label(i);
        }
    }
}

