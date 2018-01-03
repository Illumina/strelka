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

///
/// \author Chris Saunders
///

#pragma once

#include "blt_util/blt_types.hh"
#include "blt_util/pos_range.hh"
#include "starling_common/indel_core.hh"

#include <iosfwd>
#include <string>


/// key used to uniquely describe each alternate allele
///
/// pos here means zero-indexed base following the indel
/// (numerically the same as vcf POS value (+1 for 1-indexing and -1 for preceding position used in vcf)
///
struct IndelKey
{
    IndelKey(
        const pos_t p=0,
        const INDEL::index_t t=INDEL::NONE,
        const unsigned l = 0,
        const char* is = "")
        : pos(p), type(t), deletionLength(l), insertSequence(is)
    {
        validate();
    }

    // default sort is based on left-most position of the indel (note
    // we consider breakpoints to have the same left and right
    // locations)
    //
    bool
    operator<(const IndelKey& rhs) const
    {
        if (pos < rhs.pos) return true;
        if (pos != rhs.pos) return false;
        return greaterThanCore(rhs);
    }

    bool
    greaterThanCore(const IndelKey& rhs) const
    {
        if (type < rhs.type) return true;
        if (type != rhs.type) return false;
        if ((type == INDEL::NONE) ||
            (type == INDEL::BP_LEFT) ||
            (type == INDEL::BP_RIGHT)) return false;
        if (insert_length() < rhs.insert_length()) return true;
        if (insert_length() != rhs.insert_length()) return false;
        if (delete_length() < rhs.delete_length()) return true;
        if (delete_length() != rhs.delete_length()) return false;
        return (insertSequence < rhs.insertSequence);
    }

    bool
    operator==(const IndelKey& rhs) const
    {
        return ((pos == rhs.pos) &&
                (type == rhs.type) &&
                (deletionLength == rhs.deletionLength) &&
                (insertSequence == rhs.insertSequence));
    }

    pos_t right_pos() const
    {
        return pos+delete_length();
    }

    // generalize some length concepts:
    //
    unsigned
    insert_length() const
    {
        return insertSequence.size();
    }

    unsigned
    delete_length() const
    {
        return deletionLength;
    }

    const
    std::string&
    insert_seq() const
    {
        return insertSequence;
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
        if      (type == INDEL::BP_LEFT)
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

    bool
    is_breakpoint() const
    {
        return ((type == INDEL::BP_LEFT) || (type == INDEL::BP_RIGHT));
    }

    bool
    isPrimitiveDeletionAllele() const
    {
        return ((type == INDEL::INDEL) and (insertSequence.empty()) and (deletionLength>0));
    }

    bool
    isPrimitiveInsertionAllele() const
    {
        return ((type == INDEL::INDEL) and (not insertSequence.empty()) and (deletionLength==0));
    }

    bool
    isMismatch() const
    {
        return (type == INDEL::MISMATCH);
    }

    void
    validate() const;

    /// returns the global null IndelKey, any IndelKey object comparing equal to this is null
    static
    const IndelKey&
    noIndel()
    {
        static const IndelKey _noIndel;
        return _noIndel;
    }

    pos_t pos;
    INDEL::index_t type;
    unsigned deletionLength;
    std::string insertSequence; ///< insert sequence used for complete types only, not for breakends
};



#if 0
struct right_pos_indel_key_sorter
{
    bool
    operator()(
        const IndelKey& i1,
        const IndelKey& i2) const
    {
        if (i1.right_pos() < i2.right_pos()) return true;
        if (i1.right_pos() != i2.right_pos()) return false;
        return i1.greaterThanCore(i2);
    }
};
#endif


std::ostream& operator<<(std::ostream& os, const IndelKey& indelKey);
