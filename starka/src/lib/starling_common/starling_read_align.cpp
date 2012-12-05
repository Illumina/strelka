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

#include "alignment_util.hh"
#include "starling_read_align.hh"
#include "starling_read_align_clipper.hh"
#include "starling_read_align_score.hh"
#include "starling_read_align_score_indels.hh"

#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"
#include "blt_util/pos_range.hh"
#include "starling_common/indel_util.hh"

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <iostream>
#include <sstream>


//#define DEBUG_ALIGN



// Gets the lowest known min and highest known max.
//
static
known_pos_range
greatest_known_range(const known_pos_range& p1,
                     const known_pos_range& p2){
    return known_pos_range(std::min(p1.begin_pos,p2.begin_pos),
                           std::max(p1.end_pos,p2.end_pos));
}



// Check to see if an alignment overlaps any candidate indels (but not
// private indels) over the maximum range suggested by its discovery
// alignment. If at least one overlap is discovered, then the read
// goes into full re-alignment.
//
static
bool
check_for_candidate_indel_overlap(const starling_options& client_opt,
                                  const read_segment& rseg,
                                  const indel_synchronizer& isync){

#ifdef DEBUG_ALIGN
    std::cerr << "BUGBUG testing read_segment for indel overlap, sample_no: " << isync.get_sample_id() << " rseg: " << rseg;
#endif

    const pos_t seq_length(rseg.read_size());

    //get min,max bounds for each discovery alignment;
    bool is_pr_set(false);
    known_pos_range pr(0,0);
    const alignment& al(rseg.genome_align());
    if(! al.empty()) {
        pr=get_alignment_zone(al,seq_length);
        is_pr_set=true;

        // for the genomic alignment only we subtract off any edge soft-clip:
        pr.begin_pos+=apath_soft_clip_lead_size(al.path);
        pr.end_pos-=static_cast<pos_t>(apath_soft_clip_trail_size(al.path));
    }

    {
        typedef contig_align_t cat;
        const cat& ct(rseg.contig_align());
        cat::const_iterator i(ct.begin()),i_end(ct.end());
        for(;i!=i_end;++i){
            const known_pos_range k(get_alignment_zone(i->second,seq_length));
            if(! is_pr_set) {
                pr=k;
                is_pr_set=true;
            } else {
                pr=greatest_known_range(pr,k);
            }
        }
    }

    assert(is_pr_set);

#ifdef DEBUG_ALIGN
    std::cerr << "VARMIT read extends: " << pr << "\n";
#endif

    const indel_buffer& ibuff(isync.ibuff());
    const std::pair<ciiter,ciiter> ipair(ibuff.pos_range_iter(pr.begin_pos,pr.end_pos));
    for(ciiter i(ipair.first);i!=ipair.second;++i){
        const indel_key& ik(i->first);
        const indel_data& id(get_indel_data(i));
#ifdef DEBUG_ALIGN
        std::cerr << "VARMIT key: " << ik;
#endif
        // check if read intersects with indel breakpoint:
        if(! is_range_intersect_indel_breakpoints(pr,ik)) continue;
#ifdef DEBUG_ALIGN
        std::cerr << "VARMIT intersects indel: " << ik << id;
#endif

        // check if indel qualifies as candidate indel:
        if(isync.is_candidate_indel(client_opt,ik,id)){
#ifdef DEBUG_ALIGN
            std::cerr << "VARMIT read segment intersects at least one qualifying indel.\n";
#endif
            return true;
        }
    }
#ifdef DEBUG_ALIGN
    std::cerr << "VARMIT read does not intersect indel.\n";
#endif
    return false;
}



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



// transform an alignment such that any insert edge segments become
// match. insertions can be enclosed with soft-clip/hard-clip and will still be
// counted as edge insertions.
//
// segments are joined and start pos is adjusted appropriately
//
static
alignment
matchify_edge_insertion(const alignment& al) {

    using namespace ALIGNPATH;

    alignment al2;
    al2.is_fwd_strand=al.is_fwd_strand;
    al2.pos=al.pos;

    const std::pair<unsigned,unsigned> ends(get_nonclip_end_segments(al.path));
    const unsigned as(al.path.size());
    for(unsigned i(0);i<as;++i){
        const path_segment& ps(al.path[i]);
        const bool is_edge_segment((i==ends.first) || (i==ends.second));
        const bool is_clip_type(ps.type==INSERT);
        const bool is_edge_clip(is_edge_segment && is_clip_type);
        if(is_clip_type && (i==ends.first)) al2.pos-=ps.length;
        if(is_edge_clip || (ps.type==MATCH)){
            if((! al2.empty()) && (al2.path.back().type==MATCH)){
                al2.path.back().length+=ps.length;
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




// transform an alignment such that any soft-clipped edge segments
// become match.
//
// segments are joined and start pos is adjusted appropriately
//
static
alignment
matchify_edge_soft_clip(const alignment& al) {

    using namespace ALIGNPATH;

    alignment al2;
    al2.is_fwd_strand=al.is_fwd_strand;
    al2.pos=al.pos;

    const std::pair<unsigned,unsigned> ends(get_nonclip_end_segments(al.path));

    const unsigned as(al.path.size());
    for(unsigned i(0);i<as;++i){
        const path_segment& ps(al.path[i]);
        const bool is_edge_segment((i<ends.first) || (i>ends.second));
        const bool is_clip_type(ps.type==SOFT_CLIP);
        if(is_clip_type) assert(is_edge_segment);
        const bool is_edge_clip(is_edge_segment && is_clip_type);
        if(is_clip_type && (i<ends.first)) al2.pos-=ps.length;
        if(is_edge_clip || (ps.type==MATCH)){
            if((! al2.empty()) && (al2.path.back().type==MATCH)){
                al2.path.back().length+=ps.length;
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



// transform an alignment such that any soft-clipped and insert edge
// segments become match.
//
// segments are joined and start pos is adjusted appropriately
//
static
alignment
matchify_edge(const alignment& al,
              const bool ignore_soft_clip) {

    const alignment al2(matchify_edge_insertion(al));
    if(ignore_soft_clip) return al2;
    return matchify_edge_soft_clip(al2);
}



// test new alignment against the exemplar set to determine if it
// is compatible with one that already exists:
//
static
void
add_exemplar_alignment(const alignment& al,
                       const unsigned max_indel_size,
                       const bool is_genomic,
                       std::vector<alignment>& exal){

    if(! al.is_realignable(max_indel_size)) return;

    // right now we just remove soft-clipping/edge-indels and force it
    // to match
    //
    // TODO -- find a better way to handle this case:
    //
    const alignment* al_ptr(&al);
    alignment nsal;
    if(is_edge_clipped(al.path)) {
        nsal=matchify_edge(al,is_genomic);
        al_ptr=&nsal;
    }

    const unsigned es(exal.size());
    for(unsigned i(0);i<es;++i){
        if(check_and_adjust_exemplar(*al_ptr,exal[i])) return;
    }

    // no compatible alignment found!:
    exal.push_back(*al_ptr);
}



// reduce discovery alignments to a set of non-redundant exemplars:
//
static
void
get_exemplar_alignments(const read_segment& rseg,
                        const unsigned max_indel_size,
                        std::vector<alignment>& exal) {

    exal.clear();
    const alignment& al(rseg.genome_align());
    if(! al.empty()) add_exemplar_alignment(al,max_indel_size,true,exal);

    typedef contig_align_t cat;
    const cat& ct(rseg.contig_align());
    cat::const_iterator i(ct.begin()),i_end(ct.end());
    for(;i!=i_end;++i) add_exemplar_alignment(i->second,max_indel_size,false,exal);
}



static
void
dump_indel_status_map(const indel_status_map_t& ismap,
                      std::ostream& os){

    indel_status_map_t::const_iterator i(ismap.begin()), i_end(ismap.end());
    for(;i!=i_end;++i) {
        os << i->first << "status: " << i->second << "\n";
    }
}



// is an indel either a candidate indel or in at least one of the
// discovery alignments for this read?
//
static
bool
is_usable_indel(const indel_synchronizer& isync,
                const starling_options& opt,
                const indel_key& ik,
                const indel_data& id,
                const align_id_t read_id){

    return (isync.is_candidate_indel(opt,ik,id) ||
            (id.all_read_ids.count(read_id)>0) ||
            (id.tier2_map_read_ids.count(read_id)>0) ||
            (id.submap_read_ids.count(read_id)>0) ||
            (id.noise_read_ids.count(read_id)>0));
}



// find all indels in the indel_buffer which intersect a range (and
// meet candidacy/usability requirements)
static
void
add_indels_in_range(const starling_options& opt,
                    const align_id_t read_id,
                    const indel_synchronizer& isync,
                    const known_pos_range& pr,
                    indel_status_map_t& indel_status_map,
                    std::vector<indel_key>& indel_order){

    const indel_buffer& ibuff(isync.ibuff());
    const std::pair<ciiter,ciiter> ipair(ibuff.pos_range_iter(pr.begin_pos,pr.end_pos));
#ifdef DEBUG_ALIGN
    std::cerr << "VARMIT CHECKING INDELS IN RANGE: " << pr << "\n";
#endif
    for(ciiter i(ipair.first);i!=ipair.second;++i){
        const indel_key& ik(i->first);
        // check if read intersects with indel and indel is usable by this read:
#ifdef DEBUG_ALIGN
        std::cerr << "VARMIT INDEL CANDIDATE " << ik;
        std::cerr << "Intersect?: " << is_range_intersect_indel_breakpoints(pr,ik) << "\n";
        std::cerr << "Usable?: " <<  is_usable_indel(isync,opt,ik,get_indel_data(i),read_id) << "\n";
        std::cerr << "Count: " << indel_status_map.count(ik) << "\n";
#endif
        if(! is_range_intersect_indel_breakpoints(pr,ik)) continue;

        const indel_data& id(get_indel_data(i));
        if(is_usable_indel(isync,opt,ik,id,read_id) &&
           (indel_status_map.count(ik)==0) ){
            indel_status_map[ik] = false;
            indel_order.push_back(ik);
        }
    }
}



// construct an alignment which includes all of the indels turned on
// in the indel set, holding the start_position fixed to the target
// value -- indel sets should be pre-filtered for cases where an indel
// crosses the start pos, so this is treated as an error condition:
//
static
candidate_alignment
make_start_pos_alignment(const pos_t ref_start_pos,
                         const pos_t read_start_pos,
                         const bool is_fwd_strand,
                         const unsigned read_length,
                         const indel_set_t& indels){

    using namespace ALIGNPATH;

    assert(read_length>0);
    assert(ref_start_pos>=0);
    assert(read_start_pos>=0);

    candidate_alignment cal;
    cal.al.pos=ref_start_pos;
    cal.al.is_fwd_strand=is_fwd_strand;

    pos_t ref_head_pos(ref_start_pos);
    pos_t read_head_pos(read_start_pos);

    indel_set_t::const_iterator i(indels.begin()),i_end(indels.end());

    path_t& apath(cal.al.path);

    for(;i!=i_end;++i){
        const indel_key& ik(*i);

        // check that indel is a possible interception:
        const bool is_leading_read(read_start_pos!=0);
        const pos_t swap_pos(ik.pos+ik.swap_dlength);
        if((swap_pos < ref_start_pos) || ((swap_pos == ref_start_pos) && (! is_leading_read))) continue;

        // deal with leading indel, swap or right breakpoint:
        const bool is_leading_indel((swap_pos == ref_start_pos) && is_leading_read);
        if(is_leading_indel){
            assert((apath.size()==0) && (ref_head_pos==ref_start_pos));
            assert((ik.type == INDEL::INSERT) ||
                   (ik.type == INDEL::SWAP) ||
                   (ik.type == INDEL::BP_RIGHT));

            if((ik.type == INDEL::INSERT) ||
               (ik.type == INDEL::SWAP)) {
                assert(static_cast<pos_t>(ik.length)>=read_start_pos);
            }

            apath.push_back(path_segment(INSERT,read_start_pos));
            cal.leading_indel_key=ik;
            continue;
        }

        // no more leading insertion indels -- deal with regular case:
        assert((apath.size()!=0) || (read_start_pos==0));

        // note this relies on the single extra base of separation
        // required between indels during indel conflict detection:
        assert(ik.pos > ref_head_pos);

        const unsigned match_segment(ik.pos-ref_head_pos);

        // remaining read segment match segment added after indel loop:
        if(read_head_pos+match_segment>=read_length) break;

        apath.push_back(path_segment(MATCH,match_segment));
        ref_head_pos += match_segment;
        read_head_pos += match_segment;

        if       (ik.type==INDEL::INSERT ||
                  ik.type==INDEL::SWAP) {
            unsigned insert_length(read_length-read_head_pos);
            const bool is_final(ik.length>=insert_length);
            insert_length=std::min(ik.length,insert_length);
            apath.push_back(path_segment(INSERT,insert_length));
            read_head_pos += insert_length;
            if(is_final) {
                cal.trailing_indel_key=ik;
                break;
            } else {
                if(ik.type==INDEL::SWAP) {
                    apath.push_back(path_segment(DELETE,ik.swap_dlength));
                    ref_head_pos += ik.swap_dlength;
                }
            }
        } else if(ik.type==INDEL::DELETE) {
            apath.push_back(path_segment(DELETE,ik.length));
            ref_head_pos += ik.length;
        } else if(ik.type==INDEL::BP_LEFT) {
            const unsigned overhang_length(read_length-read_head_pos);
            apath.push_back(path_segment(INSERT,overhang_length));
            read_head_pos+=overhang_length;
            cal.trailing_indel_key=ik;
            break;
        } else {
            std::ostringstream oss;
            oss << "Unexpected indel state: " << INDEL::get_index_label(ik.type) << " at: " << __FILE__  << ":" << __LINE__ ;
            throw blt_exception(oss.str().c_str());
        }
    }

    assert(read_head_pos<=static_cast<pos_t>(read_length));
    if(read_head_pos<static_cast<pos_t>(read_length)) {
        apath.push_back(path_segment(MATCH,(read_length-read_head_pos)));
    }

    return cal;
}



// work backwards from end_pos to get start_pos and read_start_pos
// when the current indel set included, and then use the
// make_start_pos_alignment routine.
//
static
void
get_end_pin_start_pos(const indel_set_t& indels,
                      const unsigned read_length,
                      const pos_t ref_end_pos,
                      const pos_t read_end_pos,
                      pos_t& ref_start_pos,
                      pos_t& read_start_pos){

    assert(read_length>0);
    assert(ref_end_pos>0);
    assert(read_end_pos>0);

    ref_start_pos=ref_end_pos;
    read_start_pos=read_end_pos;

    bool is_first(true);

    // having trouble with normal reverse_iterator for this data
    // structure, so reversal is done by hand:
    indel_set_t::const_iterator i(indels.end()),i_begin(indels.begin());
    while(i!=i_begin){
        --i;
        const indel_key& ik(*i);

        // check that indel is a possible interception:
        const bool is_trailing_read(read_end_pos != static_cast<pos_t>(read_length));
        if((ik.pos > ref_end_pos) || ((ik.pos==ref_end_pos) && (! is_trailing_read))) continue;

        const bool is_trailing_indel((ik.pos == ref_end_pos) && is_trailing_read);
        if(is_trailing_indel){ // deal with trailing-edge insert/breakpoint case first
            assert((is_first) && (ref_start_pos==ref_end_pos));
            assert((ik.type == INDEL::INSERT) ||
                   (ik.type == INDEL::SWAP) ||
                   (ik.type == INDEL::BP_LEFT));

            if((ik.type == INDEL::INSERT) ||
               (ik.type == INDEL::SWAP)) {
                assert(ik.length>=(read_length-read_end_pos));
            }
        } else { // deal with normal case:
            if(is_first && (read_end_pos!=static_cast<pos_t>(read_length))){
                    log_os << "ERROR: is_first: " << is_first
                           << " read_end_pos: " << read_end_pos
                           << " read_length: " << read_length << "\n";
                    assert(0);
            }
            // note this relies on the single extra base of separation
            // required between indels during indel conflict detection:
            assert(ik.right_pos() < ref_start_pos);

            const unsigned match_segment(std::min((ref_start_pos - ik.right_pos()),read_start_pos));

            ref_start_pos -= match_segment;
            read_start_pos -= match_segment;

            if(read_start_pos==0) return;

            if       (ik.type==INDEL::INSERT || ik.type==INDEL::SWAP) {
                if(static_cast<pos_t>(ik.length) >= read_start_pos) return;
                read_start_pos -= ik.length;
                ref_start_pos -= ik.swap_dlength;
            } else if(ik.type==INDEL::DELETE) {
                ref_start_pos -= ik.length;
            } else if(ik.type==INDEL::BP_RIGHT) {
                return;
            } else {
                std::ostringstream oss;
                oss << "Unexpected indel state: " << INDEL::get_index_label(ik.type) << " at: " << __FILE__  << ":" << __LINE__ ;
                throw blt_exception(oss.str().c_str());
            }
        }
        is_first=false;
    }

    assert(read_start_pos >= 0);
    ref_start_pos -= read_start_pos;
    read_start_pos = 0;
}



struct mca_warnings {
    mca_warnings() : origin_skip(false), max_toggle_depth(false) {}
    bool origin_skip;
    bool max_toggle_depth;
};



// Recursively build potential alignment paths and push them into the
// candidate alignment set:
//
static
void
make_candidate_alignments(const starling_options& client_opt,
                          const starling_deriv_options& client_dopt,
                          const align_id_t read_id,
                          const unsigned read_length,
                          const indel_synchronizer& isync,
                          std::set<candidate_alignment>& cal_set,
                          mca_warnings& warn,
                          indel_status_map_t indel_status_map,
                          std::vector<indel_key> indel_order,
                          const unsigned depth,
                          const unsigned toggle_depth, // total number of changes made to the exemplar alignment
                          known_pos_range read_range,
                          int max_read_indel_toggle,
                          const candidate_alignment cal){
#ifdef DEBUG_ALIGN
    std::cerr << "VARMIT starting MCA depth: " << depth << "\n";
    std::cerr << "\twith cal: " << cal;
#endif

    // first step is to check for new indel overlaps and extend the
    // indel_status_map as necessary:
    //
    // note that we search for new indels in the range expanded from
    // the last alignment, but include an additional base from the
    // previous range so that we correctly overlap all potential new
    // indels.
    //
    bool is_new_indels(toggle_depth==0);
    {
        const unsigned start_ism_size(indel_status_map.size());
        const known_pos_range pr(get_soft_clip_alignment_range(cal.al));
        if(pr.begin_pos < read_range.begin_pos){
            add_indels_in_range(client_opt,read_id,isync,
                                known_pos_range(pr.begin_pos,read_range.begin_pos+1),
                                indel_status_map,indel_order);
            read_range.begin_pos = pr.begin_pos;
        }
        if(pr.end_pos > read_range.end_pos) {
            add_indels_in_range(client_opt,read_id,isync,
                                known_pos_range(read_range.end_pos-1,pr.end_pos),
                                indel_status_map,indel_order);
            read_range.end_pos = pr.end_pos;
        }

        if(! is_new_indels) {
            is_new_indels=(start_ism_size!=indel_status_map.size());
        }
    }

    // next check for recursive termination:
    if(depth == indel_order.size()){
        cal_set.insert(cal);
        return;
    }

    // each recursive step invokes (up to) 3 paths:
    //  1) is the current state of the active indel
    //  2) is the alternate state of the active indel with the start position pinned
    //  3) is the alternate state of the active indel with the end position pinned
    //
    // toggle will not be invoked for unusable indels (but those are
    // already filtered out of the list)
    //
    // ??? -- no longer true??? -- toggle will not be invoked if they
    // lead to indel conflicts with a previous indel
    //
    // start or end position may be skipped if a deletion spans one or
    // both of these points
    //
    // edge-indels only be pinned on one side
    //
    const indel_key& cindel(indel_order[depth]);

    const bool is_cindel_on(indel_status_map[cindel]);

    // alignment 1) --> unchanged case:
    make_candidate_alignments(client_opt,client_dopt,read_id,read_length,isync,cal_set,warn,
                              indel_status_map,indel_order,depth+1,toggle_depth,read_range,
                              max_read_indel_toggle,cal);

    // check whether this indel would interfere with an indel that's
    // already been flipped on:
    //
    if(! is_cindel_on){
        for(unsigned i(0);i<depth;++i){
            const indel_key& ik(indel_order[i]);
            if(indel_status_map[ik] && is_indel_conflict(ik,cindel)) return;
        }
    }

    if(is_new_indels) {
        // Check for very high candidate indel density. If found, indel
        // search toggling is turned down to the minimum level which still
        // allows simple calls (distance 1 from exemplar). The intention
        // is to allow basic indel calling to procede in practical time
        // through regions with very high numbers of candidate indels.
        //
        // This filter was introduced to handle a large set of predictions
        // made by grouper in a low complexity CT-rich region of human
        // NCBI37:chr16:33954000-33956000
        //
        // TODO: Note the expression for indel density doesn't account for
        // a possibly large number of indels intersected by a large
        // deletion in a read. This works well enough for now, but if the
        // max indel size is ever run around the order of 10k or more this
        // might start to spuriously engage the filter.
        //
        const double max_indels(read_length*client_opt.max_candidate_indel_density);
        if(indel_status_map.size()>max_indels){
            max_read_indel_toggle=1;
        } else {
            max_read_indel_toggle=client_opt.max_read_indel_toggle;
        }

        // a new stronger complexity limit:
        //
        {
            const int max_toggle(client_dopt.sal.get_max_toggle(indel_status_map.size()));
            max_read_indel_toggle=std::min(max_read_indel_toggle,max_toggle);
        }
    }

    // check whether toggling this indel would exceed the maximum
    // number of toggles made to the exemplar alignment (this is only
    // here to prevent a combinatorial blowup)
    //
    if(static_cast<int>(toggle_depth+1)>max_read_indel_toggle){
        warn.max_toggle_depth=true;
        return;
    }

    // changed cases:
    indel_status_map[cindel]=(! is_cindel_on);

    // extract only those indels that are present in the next
    // alignment:
    //
    indel_set_t current_indels;
    {
        typedef indel_status_map_t::const_iterator siter;
        const siter i_begin(indel_status_map.begin());
        const siter i_end(indel_status_map.end());

        for(siter i(i_begin);i!=i_end;++i){ if(i->second) current_indels.insert(i->first); }
    }
    // a pin on either end of the alignment is not possible/sensible
    // if:
    //
    // A) a deletion is being added which spans the pin site
    // B) an edge insertion/breakpoint is being removed from the pinned side

    // alignment 2) -- insert or delete indel and pin the start position
    //
    // test for conditions where the start pin is not possible:
    //
    {
        const pos_t ref_start_pos(cal.al.pos);

        const bool is_start_pos_delete_span(cindel.open_pos_range().is_pos_intersect(ref_start_pos));
        const bool is_start_pos_insert_span(is_cindel_on && (cindel == cal.leading_indel_key));
        const bool is_start_pin_valid(! (is_start_pos_delete_span || is_start_pos_insert_span));

        if(is_start_pin_valid) {
            const pos_t read_start_pos(apath_read_lead_size(cal.al.path));
            make_candidate_alignments(client_opt,client_dopt,read_id,read_length,isync,cal_set,warn,
                                      indel_status_map,indel_order,depth+1,toggle_depth+1,read_range,
                                      max_read_indel_toggle,
                                      make_start_pos_alignment(ref_start_pos,
                                                               read_start_pos,
                                                               cal.al.is_fwd_strand,
                                                               read_length,
                                                               current_indels));
        }
    }

    // check whether this is an equal-length swap, in which case alignment 3 is unnecessary:
    if((cindel.type==INDEL::SWAP) && (cindel.length==cindel.swap_dlength)) return;

    // alignment 3) -- insert or delete indel and pin the end position
    //
    // test for conditions where end-pin is not possible:
    //
    {
        const pos_t ref_end_pos(cal.al.pos+apath_ref_length(cal.al.path));

        const bool is_end_pos_delete_span(cindel.open_pos_range().is_pos_intersect(ref_end_pos-1));
        const bool is_end_pos_insert_span(is_cindel_on && (cindel == cal.trailing_indel_key));
        const bool is_end_pin_valid(! (is_end_pos_delete_span || is_end_pos_insert_span));

        if(is_end_pin_valid) {
            // work backwards from end_pos to get start_pos and
            // read_start_pos when the current indel set included,
            // and then used the make_start_pos_alignment routine.
            const pos_t read_end_pos(read_length-apath_read_trail_size(cal.al.path));
            pos_t ref_start_pos(0);
            pos_t read_start_pos(0);
            get_end_pin_start_pos(current_indels,read_length,
                                  ref_end_pos,read_end_pos,
                                  ref_start_pos,read_start_pos);

            // gaurd against low-frequency circular chromosome event:
            if(ref_start_pos<0) {
                warn.origin_skip=true;
            } else {
                make_candidate_alignments(client_opt,client_dopt,read_id,read_length,isync,cal_set,warn,
                                          indel_status_map,indel_order,depth+1,toggle_depth+1,read_range,
                                          max_read_indel_toggle,
                                          make_start_pos_alignment(ref_start_pos,
                                                                   read_start_pos,
                                                                   cal.al.is_fwd_strand,
                                                                   read_length,
                                                                   current_indels));
            }
        }
    }
}



struct extra_path_info {
    extra_path_info() : indel_count(0), del_size(0), ins_size(0), sum_pos(0) {}

    unsigned indel_count;
    unsigned del_size;
    unsigned ins_size;
    unsigned sum_pos;
};



static
extra_path_info
get_extra_path_info(const ALIGNPATH::path_t& p){
    using namespace ALIGNPATH;

    unsigned read_pos(0);

    extra_path_info epi;
    const unsigned ps(p.size());
    for(unsigned i(0);i<ps;++i){
        const path_segment& ps(p[i]);
        if(ps.type != MATCH) epi.indel_count++;
        if(ps.type == DELETE) {
            epi.del_size += ps.length;
            epi.sum_pos += read_pos;
        }
        if(ps.type == INSERT) {
            epi.ins_size += ps.length;
            epi.sum_pos += read_pos;
        }

        if(is_segment_type_read_length(ps.type)) read_pos += ps.length;
    }

    return epi;
}



#if 0
// some rather involved tie breaking for alignments with the same score:
//
// 1. choose the alignment with the fewest indels
// 2. choose the alignment with the smallest deletion length
// 3. choose the alignment with the smallest insertion length
//
static
bool
is_first_path_preferred(const ALIGNPATH::path_t& p1,
                       const ALIGNPATH::path_t& p2) {

    const extra_path_info epi1(get_extra_path_info(p1));
    const extra_path_info epi2(get_extra_path_info(p2));

    if(epi2.indel_count < epi1.indel_count) return false;
    if(epi2.indel_count == epi1.indel_count) {
        if(epi2.del_size < epi1.del_size) return false;
        if(epi2.del_size == epi1.del_size) {
            if(epi2.ins_size < epi1.ins_size) return false;
        }
    }
    return true;
}
#endif



static
unsigned
get_candidate_indel_count(const starling_options& client_opt,
                          const indel_synchronizer& isync,
                          const candidate_alignment& cal){

    unsigned val(0);

    indel_set_t is;
    get_alignment_indels(cal,client_opt.max_indel_size,is);

    typedef indel_set_t::const_iterator siter;
    siter i(is.begin()), i_end(is.end());
    for(;i!=i_end;++i) {
        const indel_key& ik(*i);
        if(isync.is_candidate_indel(client_opt,ik)) val++;
    }
    return val;
}



// some rather involved tie breaking for alignments with the same score:
//
// 1. choose the alignment with the fewest indels
// 2. choose the alignment with the fewest private indel
// 2. choose the alignment with the smallest insertion length
// 3. choose the alignment with the smallest deletion length
//
static
bool
is_first_cal_preferred(const starling_options& client_opt,
                       const indel_synchronizer& isync,
                       const candidate_alignment& c1,
                       const candidate_alignment& c2) {

    const extra_path_info epi1(get_extra_path_info(c1.al.path));
    const extra_path_info epi2(get_extra_path_info(c2.al.path));

    if(epi2.indel_count < epi1.indel_count) return false;
    if(epi2.indel_count == epi1.indel_count) {
        const unsigned cic1(get_candidate_indel_count(client_opt,isync,c1));
        const unsigned cic2(get_candidate_indel_count(client_opt,isync,c2));
        if(cic2 > cic1) return false;
        if(cic2 == cic1) {
            if(epi2.ins_size < epi1.ins_size) return false;
            if(epi2.ins_size == epi1.ins_size) {
                if(epi2.del_size < epi1.del_size) return false;
                if(epi2.del_size == epi1.del_size) {
                    if(epi2.sum_pos < epi1.sum_pos) return false;
                }
            }
        }
    }
    return true;
}



#if 0
// make sure the cal pool contains at least one candidate indel:
static
bool
is_cal_pool_contains_candidate(const starling_options& client_opt,
                               const depth_buffer& db,
                               const indel_buffer& ibuff,
                               const cal_pool_t& max_cal_pool) {

    const unsigned n_cal(max_cal_pool.size());
    for(unsigned i(0);i<n_cal;++i){
        if(0 < get_candidate_indel_count(client_opt,db,ibuff,*(max_cal_pool[i]))) return true;
    }
    return false;
}
#endif



static
void
finish_realignment(const starling_options& client_opt,
                   read_segment& rseg,
                   const cal_pool_t& cal_pool,
                   const double path_lnp,
                   const candidate_alignment* cal_ptr){

    if(client_opt.is_clip_ambiguous_path &&
       (cal_pool.size() > 1)){
        // soft-clip off any ambiguous regions from the alignment:
        // NOTE this can result in an empty alignment!!!
        //
        const unsigned n_cal(cal_pool.size());
        unsigned best_cal_id(n_cal);
        for(unsigned i(0);i<n_cal;++i){
            if(cal_pool[i]==cal_ptr) {
                best_cal_id=i;
                break;
            }
        }
        assert(best_cal_id != n_cal);
        get_clipped_alignment_from_cal_pool(cal_pool,best_cal_id,rseg.realignment);
        if(rseg.realignment.empty()) {
            rseg.realignment=cal_ptr->al;
#ifdef DEBUG_ALIGN
            log_os << "VARMIT clipping failed -- revert to: " << rseg.realignment;
        } else {
            log_os << "VARMIT clipped alignment: " << rseg.realignment;
#endif
        }
    } else {
        rseg.realignment=cal_ptr->al;
    }
    //rseg.realign_path_lnp=static_cast<float>(path_lnp);
}



// score the candidate alignments, find the most likely and the best for realignment:
//
static
void
score_candidate_alignments(const starling_options& client_opt,
                           const reference_contig_segment& ref,
                           read_segment& rseg,
                           indel_synchronizer& isync,
                           const std::set<candidate_alignment>& cal_set,
                           std::vector<double>& cal_set_path_lnp,
                           double& max_path_lnp,
                           const candidate_alignment*& max_cal_ptr){

    // the smooth optimum alignment and alignment pool are actually
    // used for realignment, whereas the strict max_path path
    // alingmnet is reported back to the indel_scoring routines.
    //
    // the smooth set may be different from the max pool set for two
    // reasons: (1) the path score is within smooth range of the
    // optimum (2) the smooth set is restrained by one or both edges of
    // an alignment being 'pinned' for exon alignments.
    //

    const indel_buffer& ibuff(isync.ibuff());

    // pins are used to prevent exon/intron boundaries from being moved
    // during exon realignment:
    //
    const std::pair<bool,bool> edge_pin(rseg.get_segment_edge_pin());
    const bool is_pinned(edge_pin.first || edge_pin.second);

    typedef std::set<candidate_alignment>::const_iterator citer;
    const citer i_begin(cal_set.begin()),i_end(cal_set.end());
    for(citer i(i_begin);i!=i_end;++i){
        const candidate_alignment& ical(*i);
        const double path_lnp(score_candidate_alignment(client_opt,ibuff,rseg,ical,ref));

        cal_set_path_lnp.push_back(path_lnp);

#ifdef DEBUG_ALIGN
        std::cerr << "VARMIT CANDIDATE ALIGNMENT " << ical;
        std::cerr << "score: " << path_lnp << "\n";
#endif

        if(NULL!=max_cal_ptr) {
            if(path_lnp<max_path_lnp) continue;

            // TODO -- cleaner test of float equivilence (the
            // present test should be legit given the way the
            // score is calculated, but it's still not preferred)
            //
            if((path_lnp<=max_path_lnp) &&
               is_first_cal_preferred(client_opt,isync,
                                      *max_cal_ptr,ical)) continue;
        }
        max_path_lnp=path_lnp;
        max_cal_ptr=&ical;
    }

    // if there's a pin on this segment, then find out which
    // alignments are allowed and what the max path_lnp of this
    // subset is:
    //
    double max_allowed_path_lnp(max_path_lnp);
    std::vector<bool> is_cal_allowed;
    if(is_pinned) {
        const candidate_alignment* max_allowed_cal_ptr(NULL);
        is_cal_allowed.resize(cal_set.size(),true);
        const known_pos_range gen_pr(get_strict_alignment_range(rseg.genome_align()));

        unsigned n(0);
        for(citer i(i_begin);i!=i_end;++i,++n){
            const known_pos_range cal_pr(get_strict_alignment_range(i->al));
            if       (edge_pin.first && (cal_pr.begin_pos != gen_pr.begin_pos)){
                is_cal_allowed[n] = false;
            } else if(edge_pin.second && (cal_pr.end_pos != gen_pr.end_pos)){
                is_cal_allowed[n] = false;
            }
            if(! is_cal_allowed[n]) continue;

            if(NULL!=max_allowed_cal_ptr) {
                if(cal_set_path_lnp[n]<max_allowed_path_lnp) continue;
            }
            max_allowed_path_lnp=cal_set_path_lnp[n];
            max_allowed_cal_ptr=&(*i);
        }
        assert(NULL != max_allowed_cal_ptr);
    }

    // go back through the the path_lnp values and produce a pool
    // that:
    //
    // (1) possibly allows sub maximal values (if
    // is_smoothed_alignments is set)
    //
    // (2) accounts for end pins (used on exon splice sites for
    // instance)
    //
    double allowed_lnp_range(0);
    if(client_opt.is_smoothed_alignments) {
        allowed_lnp_range=client_opt.smoothed_lnp_range;
    }

    double smooth_path_lnp(0);
    const candidate_alignment* smooth_cal_ptr(NULL);
    cal_pool_t smooth_cal_pool; // set of alignments within smooth-range of max path score

    unsigned n(0);
    for(citer i(i_begin);i!=i_end;++i,++n){
        if((cal_set_path_lnp[n]+allowed_lnp_range) < max_allowed_path_lnp) continue;
        if(is_pinned && (! is_cal_allowed[n])) continue;
        const candidate_alignment& ical(*i);
        smooth_cal_pool.push_back(&ical);
        if((NULL==smooth_cal_ptr) ||
           (! is_first_cal_preferred(client_opt,isync,
                                       *smooth_cal_ptr,ical))){
            smooth_path_lnp=cal_set_path_lnp[n];
            smooth_cal_ptr=&ical;
        }
    }

    assert(NULL != smooth_cal_ptr);

#ifdef DEBUG_ALIGN
    std::cerr << "BUBBY: key,max_path_lnp,max_path: " << rseg.key() << " " << max_path_lnp << " max_cal: " << *max_cal_ptr;

    if(smooth_cal_pool.size() > 1) {
        const unsigned n_cal(smooth_cal_pool.size());
        std::cerr << "BUBBY: " << n_cal << " final alignment pool:\n";
        for(unsigned i(0);i<n_cal;++i){
             std::cerr << "BUBBY: alignment " << i << "\n" << *(smooth_cal_pool[i]);
            const known_pos_range ipr(get_strict_alignment_range(smooth_cal_pool[i]->al));
            for(unsigned j(i+1);j<n_cal;++j){
                const known_pos_range jpr(get_strict_alignment_range(smooth_cal_pool[j]->al));
                if(ipr.begin_pos==jpr.begin_pos && ipr.end_pos==jpr.end_pos) std::cerr << "COWSLIP\n";
            }
        }
    }

    if(client_opt.is_smoothed_alignments) {
        std::cerr << "BUBBY: smooth_path_lnp,smooth_path: " << smooth_path_lnp << " smooth_cal: " << *smooth_cal_ptr;
    }
#endif

    // record "best" alignment for the purpose of re-alignment --
    // optionally clip best alignment here based on multiple best
    // candidates:
    rseg.is_realigned=true;
    finish_realignment(client_opt,rseg,smooth_cal_pool,smooth_path_lnp,smooth_cal_ptr);
}





// find the most likely alignment and most likely alignment for
// each indel state for every indel in indel_status_map
//
static
void
score_candidate_alignments_and_indels(const starling_options& opt,
                                      const starling_deriv_options& dopt,
                                      const starling_sample_options& sample_opt,
                                      const reference_contig_segment& ref,
                                      read_segment& rseg,
                                      indel_synchronizer& isync,
                                      std::set<candidate_alignment>& cal_set,
                                      const bool is_incomplete_search){

    assert(! cal_set.empty());

    // (1) score all alignments in cal_set and find the max scoring
    // alignment path:
    //
    // (1b) possibly also find a "best" scoring path for re-alignment
    // which is not (quite) the max
    //
    std::vector<double> cal_set_path_lnp;
    double max_path_lnp(0);
    const candidate_alignment* max_cal_ptr(NULL);

    score_candidate_alignments(opt,ref,rseg,isync,cal_set,cal_set_path_lnp,max_path_lnp,max_cal_ptr);


    //
    // Realignment for snp-calling and visualization is complete here,
    // remaining task is to evaluate alternate alignments as required
    // by the indel calling model.
    //

    // First determine if we even need to continue
    //

    if(! opt.is_call_indels()) return;

    // if calling indels, we only need reads with tier1 and tier2
    // mappings. if tier2 mappings aren't being used, they won't be in
    // the data, so there's no reason to check for usage here:
    //
    if(! rseg.is_treated_as_anytier_mapping()) return;

    score_indels(opt,dopt,sample_opt,rseg,isync,cal_set,is_incomplete_search,cal_set_path_lnp,max_path_lnp,max_cal_ptr);
}



static
void
get_exemplar_candidate_alignments(const starling_options& opt,
                                  const starling_deriv_options& dopt,
                                  const read_segment& rseg,
                                  const indel_synchronizer& isync,
                                  const alignment& exemplar,
                                  mca_warnings& warn,
                                  std::set<candidate_alignment>& cal_set){

    const unsigned read_length(rseg.read_size());

    indel_status_map_t indel_status_map;
    std::vector<indel_key> indel_order;

    candidate_alignment cal;
    cal.al=exemplar;

#ifdef DEBUG_ALIGN
    std::cerr << "VARMIT starting search from exemplar: " << cal;
#endif

    // Get indel set and indel order for the exemplar alignment:
    const known_pos_range pr(get_strict_alignment_range(cal.al));
    add_indels_in_range(opt,rseg.id(),isync,pr,indel_status_map,indel_order);

#ifdef DEBUG_ALIGN
    std::cerr << "VARMIT exemplar alignment range: " << pr << "\n";
#endif

    // Mark the indels which are already included in the discovery
    // alignment. Note this won't work for edge indels but we've
    // removed those from the exemplars already.
    //
    {
        indel_set_t cal_indels;
        get_alignment_indels(cal,opt.max_indel_size,cal_indels);

        typedef indel_set_t::const_iterator siter;
        const siter i_begin(cal_indels.begin());
        const siter i_end(cal_indels.end());
        for(siter i(i_begin);i!=i_end;++i) {
            if(indel_status_map.find(*i)==indel_status_map.end()) {
                log_os << "ERROR: Exemplar alignment contains indel not found in the overlap indel set\n"
                       << "\tIndel: " << *i
                       << "Exemplar overlap set:\n";
                dump_indel_status_map(indel_status_map,log_os);
                exit(EXIT_FAILURE);
            }
            indel_status_map[*i] = true;
        }
    }

    // to prevent incompatible alignments, we must put all active indels first in the order list:
    //
    typedef indel_status_map_t::const_iterator siter;
    const siter i_begin(indel_status_map.begin());
    const siter i_end(indel_status_map.end());

    indel_order.clear();
    for(siter i(i_begin);i!=i_end;++i){ if(i->second) indel_order.push_back(i->first); }
    for(siter i(i_begin);i!=i_end;++i){ if(! i->second) indel_order.push_back(i->first); }

#ifdef DEBUG_ALIGN
    std::cerr << "VARMIT exemplar starting indel_status_map:\n";
    dump_indel_status_map(indel_status_map,std::cerr);

    {
        std::cerr << "VARMIT exemplar starting indel_order:\n";
        const unsigned foo(indel_order.size());
        for(unsigned j(0);j<foo;++j) {
            std::cerr << "no: " << j << " " << indel_order[j] << "\n";
        }
    }
#endif

    // to handle soft-clip and hard-clip in the genomic alignment, we take the
    // soft-clip portion off of the alignment before entering
    // make_candidate_alignments, and add it back in afterwards.
    //
    // TODO -- something less hacky to handle soft-clip
    //
    unsigned cal_read_length(read_length);
    unsigned hc_lead(0);
    unsigned hc_trail(0);
    unsigned sc_lead(0);
    unsigned sc_trail(0);
    const bool is_exemplar_clip(is_clipped(cal.al.path));
    if(is_exemplar_clip){
        // exemplar clip condition should only be true for the
        // genomic alignment, which comes first in the exemplar list:
        //
        assert(cal_set.empty());

        apath_clip_clipper(cal.al.path,
                           hc_lead,
                           hc_trail,
                           sc_lead,
                           sc_trail);

        assert(cal_read_length >= (sc_lead+sc_trail));
        cal_read_length-=(sc_lead+sc_trail);
    }

    // launch recursive re-alignment routine starting from the current exemplar alignment:
    static const unsigned start_depth(0);
    static const unsigned start_toggle_depth(0);
    make_candidate_alignments(opt,dopt,rseg.id(),cal_read_length,isync,cal_set,warn,
                              indel_status_map,indel_order,start_depth,start_toggle_depth,pr,
                              opt.max_read_indel_toggle,cal);

    if(is_exemplar_clip) {
        // un soft-clip candidate alignments:
        std::set<candidate_alignment> cal_set2(cal_set);
        cal_set.clear();
        typedef std::set<candidate_alignment>::iterator xiter;
        const xiter i_begin(cal_set2.begin()),i_end(cal_set2.end());
        for(xiter i(i_begin);i!=i_end;++i){
            candidate_alignment cal(*i);
            apath_clip_adder(cal.al.path,
                             hc_lead,
                             hc_trail,
                             sc_lead,
                             sc_trail);
            cal_set.insert(cal);
        }
    }
}



// search for optimal realignments of the read and score alternate
// indel states in preperation for indel genotype calling
//
void
realign_and_score_read(const starling_options& opt,
                       const starling_deriv_options& dopt,
                       const starling_sample_options& sample_opt,
                       const reference_contig_segment& ref,
                       read_segment& rseg,
                       indel_synchronizer& isync){

    if(! rseg.is_valid()) {
        log_os << "ERROR: invalid alignment path associated with read segment:\n" << rseg;
        exit(EXIT_FAILURE);
    }

    // check that there are any candidate indels within bounds of the
    // discovery alignments for this read:
    //
    if(! check_for_candidate_indel_overlap(opt,rseg,isync)) return;


    // Reduce discovery alignments to a set of non-redundant exemplars.
    //
    // For the exemplar set, all contig alignments have soft-clipped
    // ends forced to match. For genomic alignments soft-clip is
    // preserved on the edge of a read, but not edge insertions.
    //
    // TODO: further adapt this to allow for actual soft-clip cases
    // (examples: read origin crossover, end of an exon... etc), by
    // finding some way for contig alignments to preserve soft-clip.
    //
    std::vector<alignment> exemplars;
    get_exemplar_alignments(rseg,opt.max_indel_size,exemplars);

    if(exemplars.empty()) return;

    // If any exemplars map to a negative start position, then skip realignment
    //
    // Note that even though such alignments will not occur in BAM, they
    // can still be produced by grouper.
    //
    {
        const unsigned n_ex(exemplars.size());
        for(unsigned i(0);i<n_ex;++i){
            if(exemplars[i].pos<0) return;
        }
    }

    // run recursive alignment search starting from each discovery
    // alignment, produce a list of candidate alignments from this
    // search
    //
    // scheme:
    //
    // 1) starting indel set are just those overlapped by the
    // discovery alignment
    //
    // 2) order is defined over the indel set -- indels present in the
    // exemplar must come first to prevent interference under this
    // scheme
    //
    // 3) when a new candidate alignment overlaps a new indel (ie. not
    // in the indel set), that indel is added to the indel set and is
    // pushed onto the end of the indel order.
    //
    // 4) indel status is recorded in the indel_status_map -- toggled
    // on or off to indicate the current state
    //
    // 5) alignments are recorded when recursion depth == indel set
    // size -- note that for this reason the indel off state has to
    // interpreted as "indel not present" rather than "reference", so
    // that all indels can be visited even if some conflict.
    //
    std::set<candidate_alignment> cal_set;
    mca_warnings warn;

    const unsigned n_ex(exemplars.size());
    for(unsigned i(0);i<n_ex;++i){
        get_exemplar_candidate_alignments(opt,dopt,rseg,isync,exemplars[i],warn,cal_set);
    }

    assert(! cal_set.empty());

    const bool is_incomplete_search(warn.origin_skip || warn.max_toggle_depth);

    // the max_toggle event is too common in genomic resequencing to have a
    // default warning:
    //
    const bool is_max_toggle_warn_enabled(opt.verbosity >= LOG_LEVEL::ALLWARN);

    if(warn.origin_skip || (warn.max_toggle_depth && is_max_toggle_warn_enabled)){
        static const char* osr="alignments crossed chromosome origin";
        static const char* mir="exceeded max number of indel switches";
        const char* reason(warn.origin_skip ? osr : mir );


        log_os << "WARNING: re-alignment skipped some alternate alignments for read: "  << rseg.key()
            //               << " at buffer_pos: " << (sread.buffer_pos+1) << "\n"
               << "\treason: " << reason << "\n";
    }

    score_candidate_alignments_and_indels(opt,dopt,sample_opt,ref,rseg,isync,cal_set,is_incomplete_search);
}
