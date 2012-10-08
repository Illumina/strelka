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

#ifndef __POS_RANGE_HH
#define __POS_RANGE_HH

#include "blt_util/blt_types.hh"

#include <algorithm>
#include <iosfwd>

///
/// utility to handle (potentially open) number ranges and test
/// intersections to either positions or other ranges.
///
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the begining of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///
/// any non-range pos value is assumed to be zero-indexed
///

struct pos_range {
    
    pos_range() : is_begin_pos(false), is_end_pos(false), begin_pos(0), end_pos(0) {}

    pos_range(const pos_t bp,const pos_t ep) 
        :  is_begin_pos(true), is_end_pos(true), begin_pos(bp), end_pos(ep) {}

    void
    clear() {
        is_begin_pos=false;
        is_end_pos=false;
        begin_pos=0;
        end_pos=0;
    }

    void
    set_begin_pos(const pos_t pos) {
        begin_pos=pos;
        is_begin_pos=true;
    }

    void
    set_end_pos(const pos_t pos) {
        end_pos=pos;
        is_end_pos=true;
    }

    bool
    is_empty() const {
        return ! (is_begin_pos || is_end_pos);
    }

    bool
    is_complete() const {
        return (is_begin_pos && is_end_pos);
    }

    inline
    bool
    is_pos_intersect(const pos_t pos) const {

        return (((! is_begin_pos) || (pos >= begin_pos)) &&
                ((! is_end_pos) || (pos < end_pos)));
    }

    bool
    is_range_intersect(const pos_range& pr) const {
        return (((! pr.is_end_pos) || (! is_begin_pos) || (pr.end_pos > begin_pos)) && 
                ((! pr.is_begin_pos) || (! is_end_pos) || (pr.begin_pos < end_pos)));
    }

    /// does this range completely overlap pr?
    bool
    is_superset_of(const pos_range& pr) const {
        return 
            (((! is_end_pos) || 
              ( pr.is_end_pos && (pr.end_pos <= end_pos) )) && 
             ((! is_begin_pos) ||
              ( pr.is_begin_pos && (pr.begin_pos >= begin_pos) )));
    }

    unsigned
    size() const {
        if(! is_complete()) return 0;
        return std::max(0,end_pos-begin_pos);
    }

    bool is_begin_pos;
    bool is_end_pos;
    pos_t begin_pos;
    pos_t end_pos;
};


// More restricted closed form. This object allows functions to express a
// closed interval requirement
//
struct known_pos_range : public pos_range {

    known_pos_range(const pos_t bp,const pos_t ep) : pos_range(bp,ep) {}

private:
    void clear();
};


std::ostream& operator<<(std::ostream& os, const pos_range& pr);



#endif
