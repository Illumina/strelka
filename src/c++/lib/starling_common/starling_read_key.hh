//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

/// \file
/// \author Chris Saunders
///

#pragma once

#include "htsapi/bam_record.hh"
#include "starling_common/starling_types.hh"

#include <iosfwd>


// information required to uniquely identify a read:
//
struct read_key
{
    explicit
    read_key(const bam_record& br) : _br_ptr(&br) {}

    int
    read_no() const
    {
        return _br_ptr->read_no();
    }

    const char*
    qname() const
    {
        return _br_ptr->qname();
    }

    bool
    operator<(const read_key& rhs) const
    {
        if (read_no()<rhs.read_no()) return true;
        if (read_no()!=rhs.read_no()) return false;
        return (strcmp(qname(),rhs.qname())<0);
    }

    bool
    operator==(const read_key& rhs) const
    {
        return ((read_no()==rhs.read_no()) and ((0==strcmp(qname(),rhs.qname()))));
    }

private:
    const bam_record* _br_ptr;
};


// debugging output:
std::ostream& operator<<(std::ostream& os, const read_key& rk);
