// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
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
        if (read_no()<rhs.read_no()) return true;
        if (read_no()==rhs.read_no()) {
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
