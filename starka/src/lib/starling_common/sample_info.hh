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
