// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
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

///
/// \author Chris Saunders
///

#pragma once

#include "blt_util/blt_types.hh"
#include "blt_util/pos_range.hh"
#include "starling_common/indel_core.hh"

#include <iosfwd>


// pos here means zero-indexed base following the indel (numerically the same as vcf POS value (+1 for 1-indexing and -1 for preceding position used in vcf)
//
// policy (for now) is that two indels which are the same except for
// their sorted sequence are treated as the same, the insert sequence
// is the first sequence encountered:
//
// length is an overloaded term:
//
// if type is insert it is the inserted sequence length
// if type is delete it is the deletion length
// if type is breakpoint it is the length of unaligned sequence stored for the other side of the breakpoint
// it type is swap it is the inserted sequence, in this case swapd_length is used to indicated the deletion length
//
struct IndelKey
{
    IndelKey(
        const pos_t p=0,
        const INDEL::index_t t=INDEL::NONE,
        const unsigned l=0,
        const unsigned sl=0)
        : pos(p), type(t), length(l), swap_dlength(sl)
    {}

    // default sort is based on left-most position of the indel (note
    // we consider breakpoints to have the same left and right
    // locations)
    //
    bool
    operator<(const IndelKey& rhs) const
    {
        if (pos < rhs.pos) return true;
        if (pos != rhs.pos) return false;
        return gtcore(rhs);
    }

    bool
    gtcore(const IndelKey& rhs) const
    {
        if (type < rhs.type) return true;
        if (type != rhs.type) return false;
        if ((type == INDEL::NONE) ||
            (type == INDEL::BP_LEFT) ||
            (type == INDEL::BP_RIGHT)) return false;
        if (length < rhs.length) return true;
        if (length != rhs.length) return false;
        return (swap_dlength < rhs.swap_dlength);
    }

    bool
    operator==(const IndelKey& rhs) const
    {
        return ((pos == rhs.pos) &&
                (type == rhs.type) &&
                (length == rhs.length) &&
                (swap_dlength == rhs.swap_dlength));
    }

    pos_t right_pos() const
    {
        if     (type==INDEL::DELETE)
        {
            return pos+length;
        }
        else if (type==INDEL::SWAP)
        {
            return pos+swap_dlength;
        }
        return pos;
    }


    // generalize some length concepts:
    //
    unsigned
    insert_length() const
    {
        if ((type == INDEL::INSERT) ||
            (type == INDEL::SWAP))
        {
            return length;
        }
        else
        {
            return 0;
        }
    }

    unsigned
    delete_length() const
    {
        if       (type == INDEL::DELETE)
        {
            return length;
        }
        else if (type == INDEL::SWAP)
        {
            return swap_dlength;
        }
        else
        {
            return 0;
        }
    }

    // correct pos range to use when we view sv's as breakpoints:
    known_pos_range breakpoint_pos_range() const
    {
        return known_pos_range(pos,right_pos());
    }

    // correct pos range to use when we view sv's as ranges
    // (ie. candidate indel interference within a single read:)
    pos_range open_pos_range() const
    {
        if       (type == INDEL::BP_LEFT)
        {
            pos_range pr;
            pr.set_begin_pos(pos);
            return pr;
        }
        else if (type == INDEL::BP_RIGHT)
        {
            pos_range pr;
            pr.set_end_pos(pos);
            return pr;
        }

        return breakpoint_pos_range();
    }

    bool is_breakpoint() const
    {
        return ((type == INDEL::BP_LEFT) || (type == INDEL::BP_RIGHT));
    }

    pos_t pos;
    INDEL::index_t type;
    unsigned length;
    unsigned swap_dlength;
};



#if 0
struct right_pos_indel_key_sorter
{
    bool
    operator()(const IndelKey& i1,
               const IndelKey& i2) const
    {
        if (i1.right_pos() < i2.right_pos()) return true;
        if (i1.right_pos() == i2.right_pos())
        {
            return i1.gtcore(i2);
        }
        return false;
    }
};
#endif


std::ostream& operator<<(std::ostream& os, const IndelKey& indelKey);

