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

#include "blt_util/log.hh"
#include "starling_common/alignment_util.hh"
#include "starling_common/starling_read_align_get_exemplars.hh"

#include "boost/foreach.hpp"

#include <cassert>

#include <iostream>



// return true if alignments are compatible -- the only allowed
// difference is edge soft clipping, if a match exists in the new
// alignment where the exemplar is clipped, the match will be
// transfered
//
// For now we assume both alignments have been "matchified" -- will
// need better soft-clip policy later.
//
// Without more information we have to assume edge inserts are private
// to the discovery alignment. This means that two otherwise
// compatible discovery alignments separated by an edge insertions are
// incompatible.  A next iteration of these routines could come up
// with a better private indel handling policy. Ideally we would want
// some step to predetermine what was private in the discovery
// alignment beforehand? But even non-private indels need to be
// preserved in the D.A. so that we no which indel alleles to evaluate
// for indel scoring -- given that the indel eval algorithm remains
// constrained to two alleles only.
//
static
bool
check_and_adjust_exemplar(const alignment& al1,
                          alignment& al2) {
    return (al1==al2);
}


/// convert segment_type to match if the segment exists before or after all match segments currently in the alignment
///
static
alignment
matchify_edge_segment_type(const alignment& al,
                           const ALIGNPATH::align_t segment_type,
                           const bool is_match_leading_edge = true,
                           const bool is_match_trailing_edge = true) {

    using namespace ALIGNPATH;

    assert(is_segment_type_read_length(segment_type));

    alignment al2;
    al2.is_fwd_strand=al.is_fwd_strand;
    al2.pos=al.pos;

    const std::pair<unsigned,unsigned> ends(get_match_edge_segments(al.path));
    const unsigned as(al.path.size());
    for(unsigned i(0); i<as; ++i) {
        const path_segment& ps(al.path[i]);
        const bool is_leading_edge_segment(i<ends.first);
        const bool is_trailing_edge_segment(i>ends.second);
        const bool is_target_type(ps.type==segment_type);
        const bool is_candidate_edge((is_match_leading_edge && is_leading_edge_segment) ||
                                     (is_match_trailing_edge && is_trailing_edge_segment));
        const bool is_edge_target(is_candidate_edge && is_target_type);
        if(is_edge_target && is_leading_edge_segment) al2.pos-=ps.length;
        if(is_edge_target || (ps.type==MATCH)) {
            if((! al2.path.empty()) && (al2.path.back().type == MATCH)) {
                al2.path.back().length += ps.length;
            } else {
                al2.path.push_back(ps);
                al2.path.back().type = MATCH;
            }
        } else {
            al2.path.push_back(ps);
        }
    }

    return al2;
}


// transform an alignment such that any insert edge segments become
// match. insertions can be enclosed with soft-clip/hard-clip and will still be
// counted as edge insertions.
//
// segments are joined and start pos is adjusted appropriately
//
static
alignment
matchify_edge_insertions(const alignment& al,
                         const bool is_match_leading_edge,
                         const bool is_match_trailing_edge) {

    return matchify_edge_segment_type(al, ALIGNPATH::INSERT,is_match_leading_edge,is_match_trailing_edge);
}



// transform an alignment such that any soft-clipped edge segments
// become match.
//
// segments are joined and start pos is adjusted appropriately
//
static
alignment
matchify_edge_soft_clip(const alignment& al) {

    return matchify_edge_segment_type(al, ALIGNPATH::SOFT_CLIP);
}


/// transform an alignment such that any edge indels
/// segments become match, and any edge deletions are removed
///
/// segments are joined and start pos is adjusted appropriately
///
static
alignment
matchify_edge_indel(const alignment& al,
                    const bool is_match_leading_edge,
                    const bool is_match_trailing_edge) {

    const alignment al2(remove_edge_deletions(al,is_match_leading_edge,is_match_trailing_edge));
    return matchify_edge_insertions(al2,is_match_leading_edge,is_match_trailing_edge);
}




// test new alignment against the exemplar set to determine if it
// is compatible with one that already exists:
//
static
void
add_exemplar_alignment(const alignment& al,
                       const unsigned max_indel_size,
                       const bool is_remove_leading_edge_indels,
                       const bool is_remove_trailing_edge_indels,
                       const bool is_remove_soft_clip,
                       std::vector<alignment>& exal) {

    if(! al.is_realignable(max_indel_size)) return;

    // right now we:
    // (1) optionally remove soft-clipping and force the clipped regions to match
    // (2) remove edge indels for grouper reads
    // TODO -- find a better way to handle this case:
    //
    const alignment* al_ptr(&al);
    alignment nial;
    if((is_remove_leading_edge_indels || is_remove_trailing_edge_indels) & is_edge_readref_len_segment(al.path)) {
        nial=matchify_edge_indel(*al_ptr,is_remove_leading_edge_indels,is_remove_trailing_edge_indels);
        al_ptr=&nial;
    }

    alignment nscal;
    if((is_remove_soft_clip) & is_soft_clipped(al.path)) {
        nscal=matchify_edge_soft_clip(*al_ptr);
        al_ptr=&nscal;
    }

    // check that this candidate exemplar does not already exist:
    BOOST_FOREACH(alignment& exemplar, exal) {
        if(check_and_adjust_exemplar(*al_ptr,exemplar)) return;
    }

    // no compatible alignment found! Add this alignment as a new exemplar:
    exal.push_back(*al_ptr);
}



// reduce discovery alignments to a set of non-redundant exemplars:
//
void
get_exemplar_alignments(const starling_options& opt,
                        const read_segment& rseg,
                        std::vector<alignment>& exal) {

    static const bool is_remove_assembler_soft_clip(true);
    const bool is_remove_mapper_soft_clip(opt.is_remap_input_softclip);

    exal.clear();

    // get exemplar from read mapper:
    const alignment& al(rseg.genome_align());
    if(! al.empty()) {
        const std::pair<bool,bool> end_pin(rseg.get_segment_edge_pin());
        const bool is_remove_leading_edge_indels(! end_pin.first);
        const bool is_remove_trailing_edge_indels(! end_pin.second);
        add_exemplar_alignment(al,opt.max_indel_size,
                               is_remove_leading_edge_indels,
                               is_remove_trailing_edge_indels,
                               is_remove_mapper_soft_clip,
                               exal);
    }

    // get additional exemplars from assembler:
    const contig_align_t& cat(rseg.contig_align());
    BOOST_FOREACH(const contig_align_t::value_type& idal, cat) {
        static const bool is_remove_leading_edge_indels(true);
        static const bool is_remove_trailing_edge_indels(true);
        add_exemplar_alignment(idal.second, opt.max_indel_size,
                               is_remove_leading_edge_indels,
                               is_remove_trailing_edge_indels,
                               is_remove_assembler_soft_clip,
                               exal);
    }

#ifdef DEBUG_ALIGN
    log_os << "VARMIT: Final exemplar set:\n";
    BOOST_FOREACH(const alignment& exemplar, exal) {
        log_os << "exemplar: " << exemplar;
    }
#endif
}


