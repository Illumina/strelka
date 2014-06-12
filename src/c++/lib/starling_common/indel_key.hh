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

#pragma once

#include "blt_util/blt_types.hh"
#include "blt_util/pos_range.hh"
#include "starling_common/indel_core.hh"
//#include "starling_common/starling_types.hh"

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
struct indel_key {

    indel_key(const pos_t p=0,
              const INDEL::index_t t=INDEL::NONE,
              const unsigned l=0,
              const unsigned sl=0,
              const unsigned mq=0,
              const unsigned bq=0)
        : pos(p), type(t), length(l), swap_dlength(sl),mapq_val(mq),baseq_val(bq) {}

    // default sort is based on left-most position of the indel (note
    // we consider breakpoints to have the same left and right
    // locations)
    //
    bool
    operator<(const indel_key& rhs) const {
        if (pos < rhs.pos) return true;
        if (pos == rhs.pos) {
            return gtcore(rhs);
        }
        return false;
    }

    bool
    gtcore(const indel_key& rhs) const {
        if (type < rhs.type) return true;
        if (type == rhs.type) {
            if ((type == INDEL::NONE) ||
                (type == INDEL::BP_LEFT) ||
                (type == INDEL::BP_RIGHT)) return false;
            if (length < rhs.length) return true;
            if (length == rhs.length) {
                if (swap_dlength < rhs.swap_dlength) return true;
            }
        }
        return false;
    }

    //updating data on the ranksums as the reads come in
    void
    addRanksumInfo(const int mapq, const int baseq, bool is_alt);

    bool
    operator==(const indel_key& rhs) const {
        return ((pos == rhs.pos) &&
                (type == rhs.type) &&
                (length == rhs.length) &&
                (swap_dlength == rhs.swap_dlength));
    }

    pos_t right_pos() const {
        if     (type==INDEL::DELETE) { return pos+length; }
        else if (type==INDEL::SWAP)   { return pos+swap_dlength; }
        return pos;
    }


    // generalize some length concepts:
    //
    unsigned
    insert_length() const {
        if ((type == INDEL::INSERT) ||
            (type == INDEL::SWAP)) {
            return length;
        } else {
            return 0;
        }
    }

    unsigned
    delete_length() const {
        if       (type == INDEL::DELETE) {
            return length;
        } else if (type == INDEL::SWAP) {
            return swap_dlength;
        } else {
            return 0;
        }
    }

    // correct pos range to use when we view sv's as breakpoints:
    known_pos_range breakpoint_pos_range() const {
        return known_pos_range(pos,right_pos());
    }

    // correct pos range to use when we view sv's as ranges
    // (ie. candidate indel interference within a single read:)
    pos_range open_pos_range() const {
        if       (type == INDEL::BP_LEFT) {
            pos_range pr;
            pr.set_begin_pos(pos);
            return pr;
        } else if (type == INDEL::BP_RIGHT) {
            pos_range pr;
            pr.set_end_pos(pos);
            return pr;
        }

        return breakpoint_pos_range();
    }

    bool is_breakpoint() const {
        return ((type == INDEL::BP_LEFT) || (type == INDEL::BP_RIGHT));
    }

    pos_t pos;
    INDEL::index_t type;
    unsigned length;
    unsigned swap_dlength;
    unsigned mapq_val;
    unsigned baseq_val;

};



#if 0
struct right_pos_indel_key_sorter {
    bool
    operator()(const indel_key& i1,
               const indel_key& i2) const {
        if (i1.right_pos() < i2.right_pos()) return true;
        if (i1.right_pos() == i2.right_pos()) {
            return i1.gtcore(i2);
        }
        return false;
    }
};
#endif


std::ostream& operator<<(std::ostream& os, const indel_key& ik);

