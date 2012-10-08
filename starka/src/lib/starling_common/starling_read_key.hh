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

#ifndef __STARLING_READ_KEY_HH
#define __STARLING_READ_KEY_HH

#include "blt_util/bam_record.hh"
#include "starling_common/starling_types.hh"

#include <iosfwd>


// information required to uniquely identify a read:
//
struct read_key {

    read_key(const bam_record& br) : _br_ptr(&br) {}

    int
    read_no() const { return _br_ptr->read_no(); }

    const char*
    qname() const { return _br_ptr->qname(); }

    bool
    operator<(const read_key& rhs) const {
        if(read_no()<rhs.read_no()) return true;
        if(read_no()==rhs.read_no()) {
            return (strcmp(qname(),rhs.qname())<0);
        }
        return false;
    }

    bool
    operator==(const read_key& rhs) const {
        return ((read_no()==rhs.read_no()) and ((0==strcmp(qname(),rhs.qname()))));
    }

private:
    const bam_record* _br_ptr;
};



// debugging output:
std::ostream& operator<<(std::ostream& os, const read_key& rk);


#endif
