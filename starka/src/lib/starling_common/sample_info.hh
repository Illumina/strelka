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

///
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the begining of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///

#ifndef __SAMPLE_INFO_HH
#define __SAMPLE_INFO_HH


// subclasses of this object can be passed around to specify sample
// information where required, the base-class handles the starling
// case of unlabeled single sample analysis.
//
struct sample_info {

    sample_info() {}
    virtual
    ~sample_info() {}

    virtual
    unsigned
    sample_size() const { return 1; }

    virtual
    const char* get_label(const unsigned /*i*/) const {
        return "";
    }
    
    virtual
    const char* get_prefix(const unsigned /*i*/,
                           const bool /*is_tier1*/) const {
        return "";
    }
};


#endif
