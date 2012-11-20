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

#ifndef __INDEL_HH
#define __INDEL_HH

#include "blt_util/blt_types.hh"
#include "blt_util/pos_range.hh"
#include "starling_common/starling_types.hh"

#include <iosfwd>
#include <map>
#include <string>
#include <set>
#include <vector>


//
// Breakpoints refer to both insertions and deletions which exceed
// varling's MAX_INDEL_SIZE. In this case we try to account for the
// breakpoint during snp-calling, but do not attempt to call the indel
// using the same methods used for small indels.
//
// Positions: Internally all positions are stored using zero-indexed
// position numbers.  For small indels and left breakpoints, we store
// the position of the first affected position. For right breakpoints,
// we store the first positions *after* the last affected
// position. Positions are stored in this manner so that the indels
// follow the starling range convention
//
// Large indel breakpoints are expected to be rare and to possibly
// have an unknown size (especially for insertions). Thus the indel
// length is not used in this case. To accomodate multiple different
// lengths at a single location, an advisory length may be inserted
// but it will not be used as part of any alignment calculation. For a
// breakpoint corresponding to a large deletion, "seq" is expected to
// represent the sequence on the other side of the breakpoint,
// starling will not look it up. In this way the break-point for *any*
// type of large-scale event can be accommodated for snp-calling. The
// seq convention for a breakpoint is to set the seq size equal to
// MAX_INDEL_SIZE
//
// Skips (introns) are *not* added here. They are treated as invariant
// alignment elements in the original read.
//
// The Swap type is initially used to represent combined
// insertion/deletion events, but will be used in the future for
// general alternate haplotypes.
//
// The "NONE" type is used for some indel lookup methods, because it sorts
// ahead of all other types at a given position:
//
//
namespace INDEL {
    enum index_t {
        NONE,
        INSERT,
        DELETE,
        BP_LEFT,
        BP_RIGHT,
        SWAP
    };

    inline
    const char*
    get_index_label(index_t id) {
        switch(id) {
        case NONE:        return "NONE";
        case INSERT:      return "INSERT";
        case DELETE:      return "DELETE";
        case BP_LEFT:  return "BP_LEFT";
        case BP_RIGHT: return "BP_RIGHT";
        case SWAP:        return "SWAP";
        default:          return "UNKNOWN";
        }
    }
}




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
              const unsigned sl=0)
        : pos(p), type(t), length(l), swap_dlength(sl) {}

    // default sort is based on left-most position of the indel (note
    // we consider breakpoints to have the same left and right
    // locations)
    //
    bool
    operator<(const indel_key& rhs) const {
        if(pos < rhs.pos) return true;
        if(pos == rhs.pos) {
            return gtcore(rhs);
        }
        return false;
    }

    bool
    gtcore(const indel_key& rhs) const {
        if(type < rhs.type) return true;
        if(type == rhs.type) {
            if((type == INDEL::NONE) ||
               (type == INDEL::BP_LEFT) ||
               (type == INDEL::BP_RIGHT)) return false;
            if(length < rhs.length) return true;
            if(length == rhs.length) {
                if(swap_dlength < rhs.swap_dlength) return true;
            }
        }
        return false;
    }

    bool
    operator==(const indel_key& rhs) const {
        return ((pos == rhs.pos) &&
                (type == rhs.type) &&
                (length == rhs.length) &&
                (swap_dlength == rhs.swap_dlength));
    }

    pos_t right_pos() const {
        if     (type==INDEL::DELETE) { return pos+length; }
        else if(type==INDEL::SWAP)   { return pos+swap_dlength; }
        return pos;
    }


    // generalize some length concepts:
    //
    unsigned
    insert_length() const {
        if((type == INDEL::INSERT) ||
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
        } else if(type == INDEL::SWAP) {
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
        } else if(type == INDEL::BP_RIGHT) {
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
};



#if 0
struct right_pos_indel_key_sorter {
    bool
    operator()(const indel_key& i1,
               const indel_key& i2) const {
        if(i1.right_pos() < i2.right_pos()) return true;
        if(i1.right_pos() == i2.right_pos()) {
            return i1.gtcore(i2);
        }
        return false;
    }
};
#endif



// Holds the alignment scores created by each read aligning across an
// indel which express the relative probability of the read aligning
// to the indel or the reference (or elsewhere in the genome). Used
// for indel genotyping.
//
struct read_path_scores {

    typedef float score_t;

    read_path_scores(const score_t r=0,
                     const score_t i=0,
                     const uint16_t ns=0,
                     const uint16_t rlen=0,
                     const bool is_t1=true)
        : ref(r)
        , indel(i)
        , nsite(ns)
        , read_length(rlen)
        , is_tier1_read(is_t1)
    {}

    void
    insert_alt(const indel_key& ik,
               const score_t a) {
        const unsigned ais(static_cast<unsigned>(alt_indel.size()));
        if(ais < 2) {
            alt_indel.push_back(std::make_pair(ik,a));
        } else {
            unsigned min_index(ais);
            score_t min(a);
            for(unsigned i(0);i<ais;++i) {
                if(alt_indel[i].second < min) {
                    min = alt_indel[i].second;
                    min_index = i;
                }
            }
            if(min_index<ais) {
                alt_indel[min_index] = std::make_pair(ik,a);
            }
        }
    }

    score_t ref;
    score_t indel;
    uint16_t nsite;

//    score_t alt;

    // store up to 2 highest scoring alternate indels
    typedef std::vector<std::pair<indel_key,score_t> > alt_indel_t;
    alt_indel_t alt_indel;

    // used to set expected het allele ratio:
    uint16_t read_length;

    // used to filter for/against tier2 data:
    bool is_tier1_read;
};



struct indel_data {
    indel_data()
        : is_external_candidate(false)
        , is_candidate_indel_cached(false)
    {}

    // add read and contig id evidence from another indel_data
    // structure -- this is not a copy ctor
    //
    void
    add_indel_data_evidence(const indel_data& id) {
        add_evidence(contig_ids,id.contig_ids);
        add_evidence(all_read_ids,id.all_read_ids);
        add_evidence(tier2_map_read_ids,id.tier2_map_read_ids);
        add_evidence(submap_read_ids,id.submap_read_ids);
        add_evidence(noise_read_ids,id.noise_read_ids);
        if(id.is_external_candidate) is_external_candidate=true;
    }


    std::string seq; // <- TODO: potentially this becomes a
                     // compilation of all insert sequences and some
                     // type of consensus is created, for now it's
                     // just the first sequence encountered

    // map_read_ids refers to the read(s) that support the indel
    // through their genomic alignments, given that those alignments
    // meet the snp-caller's mapping thresholds.
    //
    // submap_read_ids are genomic alignments that fall below the
    // mapping thresholds.
    //
    // all_read_ids contains the list of reads which either have a
    // genomic alignment passing the mapping criteria or have a contig
    // alignment.
    //
    // contig_ids are the grouper contig(s) which support the indel
    //
    typedef std::set<align_id_t> evidence_t;
    evidence_t contig_ids;
    evidence_t all_read_ids; // all contig and tier1 read_ids
    evidence_t tier2_map_read_ids;
    evidence_t submap_read_ids;

    // noise_read_ids indicates that an input alignment had a mismatch
    // fraction which was too high for the alignment to qualify as
    // support for the indel.
    //
    evidence_t noise_read_ids;

    // candidates can be provided from external sources as well:
    bool is_external_candidate;

    // this structure represents support for the indel among all reads
    // which cross an indel breakpoint by a sufficient margin after
    // re-alignment:
    //
    typedef std::map<align_id_t,read_path_scores> score_t;
    score_t read_path_lnp;
    evidence_t suboverlap_tier1_read_ids; // the reads which cross an indel breakpoint, but not by enough
                                          // to be entered into the scores list
    evidence_t suboverlap_tier2_read_ids;

    mutable bool is_candidate_indel_cached;
    mutable bool is_candidate_indel;

private:
    // a = a U b
    static
    void
    add_evidence(evidence_t& a,
                 const evidence_t& b) {
        a.insert(b.begin(),b.end());
    }
};


struct indel {
    indel_key key;
    indel_data data;
};


// debuging dumps:
std::ostream& operator<<(std::ostream& os, const read_path_scores& rps);

std::ostream& operator<<(std::ostream& os, const indel_key& ik);
std::ostream& operator<<(std::ostream& os, const indel_data& id);
std::ostream& operator<<(std::ostream& os, const indel& in);



#endif
