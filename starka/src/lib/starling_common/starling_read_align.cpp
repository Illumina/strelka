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

#include "alignment_util.hh"
#include "starling_read_align.hh"
#include "starling_read_align_clipper.hh"
#include "starling_read_align_get_exemplars.hh"
#include "starling_read_align_score.hh"
#include "starling_read_align_score_indels.hh"

#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"
#include "blt_util/pos_range.hh"
#include "starling_common/indel_util.hh"

#include "boost/foreach.hpp"

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <iostream>
#include <sstream>


//#define DEBUG_ALIGN

/// information associated with each candidate indel intersecting an alignment
struct starling_align_indel_info {

    starling_align_indel_info() :
        is_present(false),
        is_remove_only(false)
    {}

    bool is_present;
    bool is_remove_only; // candidate can be toggled off during search but not added
};



std::ostream&
operator<<(std::ostream& os, const starling_align_indel_info& ii) {

    os << "is_present: " << ii.is_present << " is_remove_only: " << ii.is_remove_only;
    return os;
}


typedef std::map<indel_key,starling_align_indel_info> starling_align_indel_status;


// Gets the lowest known min and highest known max.
//
static
known_pos_range
greatest_known_range(const known_pos_range& p1,
                     const known_pos_range& p2) {
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
check_for_candidate_indel_overlap(const starling_options& opt,
                                  const read_segment& rseg,
                                  const indel_synchronizer& isync) {

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

        if(! opt.is_remap_input_softclip) {
            // for the genomic alignment only we subtract off any edge soft-clip:
            pr.begin_pos+=apath_soft_clip_lead_size(al.path);
            pr.end_pos-=static_cast<pos_t>(apath_soft_clip_trail_size(al.path));
        }
    }

    {
        typedef contig_align_t cat;
        const cat& ct(rseg.contig_align());
        cat::const_iterator i(ct.begin()),i_end(ct.end());
        for(; i!=i_end; ++i) {
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
    for(ciiter i(ipair.first); i!=ipair.second; ++i) {
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
        if(isync.is_candidate_indel(opt,ik,id)) {
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


static
void
dump_indel_status(const starling_align_indel_status& ismap,
                  std::ostream& os) {

    BOOST_FOREACH(const starling_align_indel_status::value_type& is, ismap) {
        os << is.first << "status: " << is.second << "\n";
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
                const align_id_t read_id) {

    return (isync.is_candidate_indel(opt,ik,id) ||
            (id.all_read_ids.count(read_id)>0) ||
            (id.tier2_map_read_ids.count(read_id)>0) ||
            (id.submap_read_ids.count(read_id)>0) ||
            (id.noise_read_ids.count(read_id)>0));
}



/// find all indels in the indel_buffer which intersect a range (and
/// meet candidacy/usability requirements)
static
void
add_indels_in_range(const starling_options& opt,
                    const align_id_t read_id,
                    const indel_synchronizer& isync,
                    const known_pos_range& pr,
                    starling_align_indel_status& indel_status_map,
                    std::vector<indel_key>& indel_order) {

    const indel_buffer& ibuff(isync.ibuff());
    const std::pair<ciiter,ciiter> ipair(ibuff.pos_range_iter(pr.begin_pos,pr.end_pos));
#ifdef DEBUG_ALIGN
    std::cerr << "VARMIT CHECKING INDELS IN RANGE: " << pr << "\n";
#endif
    for(ciiter i(ipair.first); i!=ipair.second; ++i) {
        const indel_key& ik(i->first);
        // check if read intersects with indel and indel is usable by this read:
#ifdef DEBUG_ALIGN
        std::cerr << "VARMIT INDEL CANDIDATE " << ik;
        std::cerr << "Intersect?: " << is_range_intersect_indel_breakpoints(pr,ik) << "\n";
        std::cerr << "Usable?: " <<  is_usable_indel(isync,opt,ik,get_indel_data(i),read_id) << "\n";
        std::cerr << "Count: " << indel_status_map.count(ik) << "\n";
#endif
        if(! is_range_adjacent_indel_breakpoints(pr,ik)) continue;

        const bool is_remove_only(! is_range_intersect_indel_breakpoints(pr,ik));

        // if indel is already present, it may be possible to promote this indel from
        // adjacent to an intersection:
        if(indel_status_map.count(ik)) {
            if((! is_remove_only) && indel_status_map[ik].is_remove_only) {
                indel_status_map[ik].is_remove_only = false;
            }
        } else {
            const indel_data& id(get_indel_data(i));
            if(is_usable_indel(isync,opt,ik,id,read_id)) {
                indel_status_map[ik].is_present = false;
                indel_status_map[ik].is_remove_only = is_remove_only;
                indel_order.push_back(ik);
            }
        }
    }
}



static
void
add_path_segment(ALIGNPATH::path_t& apath,
                 ALIGNPATH::align_t inc_type,
                 pos_t& inc_pos,
                 const unsigned increment) {
    using namespace ALIGNPATH;

    apath.push_back(path_segment(inc_type,increment));
    inc_pos += increment;
}



// construct an alignment which includes all of the indels turned on
// in the indel set, holding the start_position fixed to the target
// value -- indel sets should be pre-filtered for cases where an indel
// crosses the start pos, so this is treated as an error condition:
//
/// see unit tests
static
candidate_alignment
make_start_pos_alignment(const pos_t ref_start_pos,
                         const pos_t read_start_pos,
                         const bool is_fwd_strand,
                         const unsigned read_length,
                         const indel_set_t& indels) {

    using namespace ALIGNPATH;

    assert(read_length>0);
    assert(ref_start_pos>=0);
    assert(read_start_pos>=0);

    // if true, this read contains a leading insert,swap,softclip,etc...
    const bool is_leading_read(read_start_pos!=0);

    candidate_alignment cal;
    cal.al.pos=ref_start_pos;
    cal.al.is_fwd_strand=is_fwd_strand;

    pos_t ref_head_pos(ref_start_pos);
    pos_t read_head_pos(read_start_pos);

    path_t& apath(cal.al.path);

    BOOST_FOREACH(const indel_key& ik, indels) {

        // don't consider indels which can't intersect the read:
        if(ik.right_pos() < ref_start_pos) continue;
        if((ik.right_pos() == ref_start_pos) && (! is_leading_read)) continue;

        if(apath.empty()) assert(ref_head_pos==ref_start_pos);

        // deal with leading indel, swap or right breakpoint:
        const bool is_first_intersecting_indel(apath.empty());

        if(is_leading_read && is_first_intersecting_indel) {
            if(ik.pos != ref_start_pos) {
                std::ostringstream oss;
                oss << "ERROR: anomalous condition for indel candidate: " << ik << "\n"
                    << "\tref_start_pos: " << ref_start_pos << "\n"
                    << "\tread_start_pos: " << read_start_pos << "\n"
                    << "\tref_head_pos: " << ref_head_pos << "\n"
                    << "\tread_head_pos: " << read_head_pos << "\n"
                    << "\tis_fwd_strand: " << is_fwd_strand << "\n"
                    << "\tread_length: " << read_length << "\n";

                oss << "\tfull indel set: ";
                BOOST_FOREACH(const indel_key& ik2, indels) {
                    oss << "\t\t" << ik2;
                }

                throw blt_exception(oss.str().c_str());
            }

            assert((apath.empty()) && (ref_head_pos==ref_start_pos));
            assert((ik.type == INDEL::INSERT) ||
                   (ik.type == INDEL::SWAP) ||
                   (ik.type == INDEL::BP_RIGHT));

            if((ik.type == INDEL::INSERT) ||
               (ik.type == INDEL::SWAP)) {
                assert(static_cast<pos_t>(ik.length)>=read_start_pos);
            }

            apath.push_back(path_segment(INSERT,read_start_pos));
            if(ik.type == INDEL::SWAP) {
                add_path_segment(apath,DELETE,ref_head_pos,ik.swap_dlength);
            }
            cal.leading_indel_key=ik;
            continue;
        }

        // no more leading insertion indels -- deal with regular case:
        assert((apath.size()!=0) || (read_start_pos==0));

        // note this relies on the single extra base of separation
        // required between indels during indel conflict detection:
        const bool is_edge_delete((INDEL::DELETE == ik.type) && (ik.pos == ref_start_pos));
        if((ik.pos <= ref_head_pos) && (!is_edge_delete)) {
            std::ostringstream oss;
            oss << "ERROR: indel candidate: " << ik << " is not greater than ref_head_pos: " << ref_head_pos
                << ". Cannot resolve indel with candidate read alignment: " << cal << "\n";
            throw blt_exception(oss.str().c_str());
        }

        const unsigned match_segment(ik.pos-ref_head_pos);

        assert(match_segment>0);

        // remaining read segment match segment added after indel loop:
        if(read_head_pos+match_segment>=read_length) break;

        apath.push_back(path_segment(MATCH,match_segment));
        ref_head_pos += match_segment;
        read_head_pos += match_segment;

        if       (ik.type==INDEL::INSERT ||
                  ik.type==INDEL::SWAP) {

            if(ik.type==INDEL::SWAP) {
                add_path_segment(apath,DELETE,ref_head_pos,ik.swap_dlength);
            }

            const unsigned max_insert_length(read_length-read_head_pos);
            const unsigned insert_length(std::min(ik.length,max_insert_length));
            add_path_segment(apath,INSERT,read_head_pos,insert_length);

            const bool is_final(ik.length>=max_insert_length);
            if(is_final) {
                cal.trailing_indel_key=ik;
                break;
            }
        } else if(ik.type==INDEL::DELETE) {
            add_path_segment(apath,DELETE,ref_head_pos,ik.length);
        } else if(ik.type==INDEL::BP_LEFT) {
            const unsigned overhang_length(read_length-read_head_pos);
            add_path_segment(apath,INSERT,read_head_pos,overhang_length);
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
/// see unit tests
static
void
get_end_pin_start_pos(const indel_set_t& indels,
                      const unsigned read_length,
                      const pos_t ref_end_pos,
                      const pos_t read_end_pos,
                      pos_t& ref_start_pos,
                      pos_t& read_start_pos) {

    assert(read_length>0);
    assert(ref_end_pos>0);
    assert(read_end_pos>0);

    ref_start_pos=ref_end_pos;
    read_start_pos=read_end_pos;

    // if true, read contains trailing insert, swap or open-ended event:
    const bool is_trailing_read(read_end_pos != static_cast<pos_t>(read_length));

    bool is_first(true);

    // having trouble with normal reverse_iterator for this data
    // structure, so reversal is done by hand:
    indel_set_t::const_iterator i(indels.end()),i_begin(indels.begin());
    while(i!=i_begin) {
        --i;
        const indel_key& ik(*i);

        // check that indel actually intersects the read:
        if(ik.pos > ref_end_pos) continue;
        if((ik.pos == ref_end_pos) && (! is_trailing_read)) continue;

        const bool is_trailing_indel(ik.right_pos() == ref_end_pos);

        if(is_trailing_indel) { // deal with trailing-edge insert/breakpoint case first
            assert((is_first) && (ref_start_pos==ref_end_pos));
            assert((ik.type == INDEL::INSERT) ||
                   (ik.type == INDEL::DELETE) ||
                   (ik.type == INDEL::SWAP) ||
                   (ik.type == INDEL::BP_LEFT));

            if((ik.type == INDEL::INSERT) ||
               (ik.type == INDEL::SWAP)) {
                assert(ik.length>=(read_length-read_end_pos));
            }

            if       (ik.type==INDEL::SWAP) {
                ref_start_pos -= ik.swap_dlength;
            } else if(ik.type==INDEL::DELETE) {
                ref_start_pos -= ik.length;
            }
        } else { // deal with normal case:
            if(is_first && (read_end_pos!=static_cast<pos_t>(read_length))) {
                std::ostringstream oss;
                oss << "ERROR: is_first: " << is_first
                    << " read_end_pos: " << read_end_pos
                    << " read_length: " << read_length << "\n";
                throw blt_exception(oss.str().c_str());
            }

            // note the excluding 'equals' relationship relies on the single extra base of separation
            // required between indels during indel conflict detection:
            //const bool is_edge_delete((INDEL::DELETE == ik.type) && (ik.right_pos() == ref_end_pos));

            // new indel must end at least one base below the current ref head (otherwise it would be
            // an interfering indel):
            //
            if((ik.right_pos() >= ref_start_pos)) { //&& (! is_edge_delete)) {
                std::ostringstream oss;
                oss << "Unexpected indel position: indel: " << ik;
                oss << "\tref_start_pos: " << ref_start_pos << " ref_end_pos: " << ref_end_pos << "\n";
                throw blt_exception(oss.str().c_str());
            }

            const unsigned match_segment(std::min((ref_start_pos - ik.right_pos()),read_start_pos));

            ref_start_pos -= match_segment;
            read_start_pos -= match_segment;

            if(read_start_pos==0) return;

            if       (ik.type==INDEL::INSERT || ik.type==INDEL::SWAP) {
                ref_start_pos -= ik.swap_dlength;
                if(static_cast<pos_t>(ik.length) >= read_start_pos) return;
                read_start_pos -= ik.length;
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



/// to prevent incomplete search, we must put new non-present remove_only indels at the end of the list:
static
void
sort_remove_only_indels_last(const starling_align_indel_status& indel_status_map,
                             std::vector<indel_key>& indel_order,
                             const unsigned current_depth = 0) {

    typedef std::vector<indel_key>::const_iterator siter;
    const siter i_begin(indel_order.begin());
    const siter i_new(i_begin+current_depth);
    const siter i_end(indel_order.end());

    std::vector<indel_key> indel_order2;
    for(siter i(i_begin); i!=i_new; ++i) { indel_order2.push_back(*i); }
    for(siter i(i_new); i!=i_end; ++i) {
        const starling_align_indel_info& sai(indel_status_map.find(*i)->second);
        if(  (sai.is_present || (! sai.is_remove_only))) indel_order2.push_back(*i);
    }
    for(siter i(i_new); i!=i_end; ++i) {
        const starling_align_indel_info& sai(indel_status_map.find(*i)->second);
        if(! (sai.is_present || (! sai.is_remove_only))) indel_order2.push_back(*i);
    }
    indel_order.swap(indel_order2);
}



struct mca_warnings {
    mca_warnings() : origin_skip(false), max_toggle_depth(false) {}
    bool origin_skip;
    bool max_toggle_depth;
};



static
void
add_pin_exception_info(
    const char* label,
    const unsigned depth,
    const candidate_alignment& cal,
    const candidate_alignment& start_cal,
    const pos_t ref_start_pos,
    const pos_t read_start_pos,
    const indel_key& cindel,
    const indel_set_t& current_indels)
{
    log_os << "\nException caught while building " << label << "-pinned alignment candidate at depth: " << depth << "\n"
           << "\tcal: " << cal
           << "\tstart_cal: " << start_cal
           << "\tref_start_pos: " << ref_start_pos << "\n"
           << "\tread_start_pos: " << read_start_pos << "\n"
           << "this_indel: " << cindel;
    BOOST_FOREACH(const indel_key& ik, current_indels) {
        log_os << "current_indels: " << ik;
    }
}



/// Recursively build potential alignment paths and push them into the
/// candidate alignment set:
///
static
void
make_candidate_alignments(const starling_options& client_opt,
                          const starling_deriv_options& client_dopt,
                          const align_id_t read_id,
                          const unsigned read_length,
                          const indel_synchronizer& isync,
                          std::set<candidate_alignment>& cal_set,
                          mca_warnings& warn,
                          starling_align_indel_status indel_status_map,
                          std::vector<indel_key> indel_order,
                          const unsigned depth,
                          const unsigned toggle_depth, // total number of changes made to the exemplar alignment
                          known_pos_range read_range,
                          int max_read_indel_toggle,
                          const candidate_alignment& cal) {

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
        if(pr.begin_pos < read_range.begin_pos) {
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

        if(is_new_indels) {
            sort_remove_only_indels_last(indel_status_map,indel_order,start_ism_size);
        }
    }

    // next check for recursive termination:
    if(depth == indel_order.size()) {
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
    const bool is_cindel_on(indel_status_map[cindel].is_present);

    // alignment 1) --> unchanged case:
    try {
        make_candidate_alignments(client_opt,client_dopt,read_id,read_length,isync,cal_set,warn,
                                  indel_status_map,indel_order,depth+1,toggle_depth,read_range,
                                  max_read_indel_toggle,cal);
    } catch(...) {
        log_os << "\nException caught while building default alignment candidate at depth: " << depth << "\n"
               << "\tcal: " << cal
               << "this_indel: " << cindel;
        throw;
    }

    if(! is_cindel_on) {
        // check whether this is a remove only indel:
        if(indel_status_map[cindel].is_remove_only) return;

        // check whether this indel would interfere with an indel that's
        // already been toggled on:
        //
        for(unsigned i(0); i<depth; ++i) {
            const indel_key& ik(indel_order[i]);
            if(indel_status_map[ik].is_present && is_indel_conflict(ik,cindel)) return;
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
        if(indel_status_map.size()>max_indels) {
            max_read_indel_toggle=1;
        } else {
            max_read_indel_toggle=client_opt.max_read_indel_toggle;
        }

        // a new stronger complexity limit on search based on total candidate indels crossing the read:
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
    if(static_cast<int>(toggle_depth+1)>max_read_indel_toggle) {
        warn.max_toggle_depth=true;
        return;
    }

    // changed cases:
    indel_status_map[cindel].is_present=(! is_cindel_on);

    // extract only those indels that are present in the next
    // alignment:
    //
    indel_set_t current_indels;
    BOOST_FOREACH(const starling_align_indel_status::value_type& is, indel_status_map) {
        if(is.second.is_present) current_indels.insert(is.first);
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
        const bool is_start_pos_indel_span(is_cindel_on && (cindel == cal.leading_indel_key));
        const bool is_start_pin_valid(! (is_start_pos_delete_span || is_start_pos_indel_span));

#ifdef DEBUG_ALIGN
        std::cerr << "VARMIT toggling MCA depth: " << depth << "\n";
        std::cerr << "VARMIT current indel: " << cindel;
        std::cerr << "VARMIT current indel on?: " << is_cindel_on << "\n";
        std::cerr << "VARMIT start-pin valid?: " << is_start_pin_valid << "\n";
#endif

        if(is_start_pin_valid) {
            const pos_t read_start_pos(apath_read_lead_size(cal.al.path));
            candidate_alignment start_cal;
            try {
                start_cal = make_start_pos_alignment(ref_start_pos,
                                                     read_start_pos,
                                                     cal.al.is_fwd_strand,
                                                     read_length,
                                                     current_indels);

                make_candidate_alignments(client_opt,client_dopt,read_id,read_length,isync,cal_set,warn,
                                          indel_status_map,indel_order,depth+1,toggle_depth+1,read_range,
                                          max_read_indel_toggle,start_cal);
            } catch (...) {
                add_pin_exception_info("start",depth,cal,start_cal,ref_start_pos,read_start_pos,cindel,current_indels);
                throw;
            }
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

        // end pin is not possible when
        // (1) an indel deletes through the end-pin position
        // (2) we try to remove a trailing indel [TODO seems like same rule should be in place for adding a trailing indel]
        const bool is_end_pos_delete_span(cindel.open_pos_range().is_pos_intersect(ref_end_pos-1));
        const bool is_end_pos_indel_span(is_cindel_on && (cindel == cal.trailing_indel_key));
        const bool is_end_pin_valid(! (is_end_pos_delete_span || is_end_pos_indel_span));

#ifdef DEBUG_ALIGN
        std::cerr << "VARMIT toggling MCA depth: " << depth << "\n";
        std::cerr << "VARMIT current indel: " << cindel;
        std::cerr << "VARMIT current indel on?: " << is_cindel_on << "\n";
        std::cerr << "VARMIT end-pin valid?: " << is_end_pin_valid << "\n";
#endif

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

            // guard against low-frequency circular chromosome event:
            if(ref_start_pos<0) {
                warn.origin_skip=true;
            } else {
                candidate_alignment start_cal;
                try {
                    start_cal = make_start_pos_alignment(ref_start_pos,
                                                         read_start_pos,
                                                         cal.al.is_fwd_strand,
                                                         read_length,
                                                         current_indels);

                    make_candidate_alignments(client_opt,client_dopt,read_id,read_length,isync,cal_set,warn,
                                              indel_status_map,indel_order,depth+1,toggle_depth+1,read_range,
                                              max_read_indel_toggle,start_cal);
                } catch (...) {
                    add_pin_exception_info("end",depth,cal,start_cal,ref_start_pos,read_start_pos,cindel,current_indels);
                    log_os << "ref_end_pos: " << ref_end_pos << "\n"
                           << "read_end_pos: " << read_end_pos << "\n";
                    throw;
                }
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
get_extra_path_info(const ALIGNPATH::path_t& p) {
    using namespace ALIGNPATH;

    unsigned read_pos(0);

    extra_path_info epi;
    BOOST_FOREACH(const path_segment& ps, p) {
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
                          const candidate_alignment& cal) {

    indel_set_t is;
    get_alignment_indels(cal,client_opt.max_indel_size,is);

    unsigned val(0);
    BOOST_FOREACH(const indel_key& ik, is) {
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
    for(unsigned i(0); i<n_cal; ++i) {
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
                   const double /*path_lnp*/,
                   const candidate_alignment* cal_ptr) {

    if(client_opt.is_clip_ambiguous_path &&
       (cal_pool.size() > 1)) {
        // soft-clip off any ambiguous regions from the alignment:
        // NOTE this can result in an empty alignment!!!
        //
        const unsigned n_cal(cal_pool.size());
        unsigned best_cal_id(n_cal);
        for(unsigned i(0); i<n_cal; ++i) {
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



static
bool
is_alignment_spanned_by_range(const known_pos_range pr,
                              const alignment& al)
{
    return pr.is_superset_of(get_strict_alignment_range(al));
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
                           const candidate_alignment*& max_cal_ptr)
{
    // the smooth optimum alignment and alignment pool are actually
    // used for realignment, whereas the strict max_path path
    // alignment is reported back to the indel_scoring routines.
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
    const citer cal_set_begin(cal_set.begin()), cal_set_end(cal_set.end());
    for(citer cal_iter(cal_set_begin); cal_iter!=cal_set_end; ++cal_iter) {
        const candidate_alignment& ical(*cal_iter);
        const double path_lnp(score_candidate_alignment(client_opt,ibuff,rseg,ical,ref));

        cal_set_path_lnp.push_back(path_lnp);

#ifdef DEBUG_ALIGN
        std::cerr << "VARMIT CANDIDATE ALIGNMENT " << ical;
        std::cerr << "score: " << path_lnp << "\n";
#endif

        if(NULL!=max_cal_ptr) {
            if(path_lnp<max_path_lnp) continue;

            // TODO -- cleaner test of float equivalence (the
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

        unsigned cal_index(0);
        for(citer cal_iter(cal_set_begin); cal_iter!=cal_set_end; ++cal_iter,++cal_index) {
            const known_pos_range cal_pr(get_strict_alignment_range(cal_iter->al));
            if       (edge_pin.first && (cal_pr.begin_pos != gen_pr.begin_pos)) {
                is_cal_allowed[cal_index] = false;
            } else if(edge_pin.second && (cal_pr.end_pos != gen_pr.end_pos)) {
                is_cal_allowed[cal_index] = false;
            }
            if(! is_cal_allowed[cal_index]) continue;

            if(NULL!=max_allowed_cal_ptr) {
                if(cal_set_path_lnp[cal_index]<max_allowed_path_lnp) continue;
            }
            max_allowed_path_lnp=cal_set_path_lnp[cal_index];
            max_allowed_cal_ptr=&(*cal_iter);
        }

        if(NULL == max_allowed_cal_ptr) {
            std::ostringstream oss;
            oss << "ERROR: reached anomalous state during search for most likely exon alignment.\n";
            oss << "\tread_segment: " << rseg << "\n";
            oss << "\tCandidate alignments:\n";
            for(citer cal_iter(cal_set_begin); cal_iter!=cal_set_end; ++cal_iter) {
                oss << *cal_iter << "\n";
            }
            throw blt_exception(oss.str().c_str());
        }
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

    unsigned cal_index(0);
    for(citer cal_iter(cal_set_begin); cal_iter!=cal_set_end; ++cal_iter,++cal_index) {
        if((cal_set_path_lnp[cal_index]+allowed_lnp_range) < max_allowed_path_lnp) continue;
        if(is_pinned && (! is_cal_allowed[cal_index])) continue;
        const candidate_alignment& ical(*cal_iter);
        smooth_cal_pool.push_back(&ical);
        if((NULL==smooth_cal_ptr) ||
           (! is_first_cal_preferred(client_opt,isync,
                                     *smooth_cal_ptr,ical))) {
            smooth_path_lnp=cal_set_path_lnp[cal_index];
            smooth_cal_ptr=&ical;
        }
    }

    assert(NULL != smooth_cal_ptr);

#ifdef DEBUG_ALIGN
    std::cerr << "BUBBY: key,max_path_lnp,max_path: " << rseg.key() << " " << max_path_lnp << " max_cal: " << *max_cal_ptr;

    if(smooth_cal_pool.size() > 1) {
        const unsigned n_cal(smooth_cal_pool.size());
        std::cerr << "BUBBY: " << n_cal << " final alignment pool:\n";
        for(unsigned i(0); i<n_cal; ++i) {
            std::cerr << "BUBBY: alignment " << i << "\n" << *(smooth_cal_pool[i]);
            const known_pos_range ipr(get_strict_alignment_range(smooth_cal_pool[i]->al));
            for(unsigned j(i+1); j<n_cal; ++j) {
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
                                      const bool is_incomplete_search) {

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



/// Initialize a candidate alignment with edge indels set according
/// to those in the input alignment
///
/// Normally we're worried about incomplete edge insertions.
/// In this case load them on the cal edges.
///
static
void
load_cal_with_edge_indels(const alignment& al,
                          candidate_alignment& cal) {

    using namespace ALIGNPATH;

    cal.al=al;

    pos_t ref_pos(al.pos);

    const std::pair<unsigned,unsigned> ends(get_match_edge_segments(al.path));
    const unsigned as(al.path.size());
    for(unsigned i(0); i<as; ++i) {
        const path_segment& ps(al.path[i]);
        if((INSERT == ps.type) || (DELETE == ps.type)) {
            INDEL::index_t itype;
            if(INSERT == ps.type) { itype = INDEL::INSERT; }
            else                  { itype = INDEL::DELETE; }
            const indel_key ik(ref_pos,itype,ps.length);
            if     (i<ends.first)  { cal.leading_indel_key = ik; }
            else if(i>ends.second) { cal.trailing_indel_key = ik; }
        }
        if(is_segment_type_ref_length(ps.type)) ref_pos += ps.length;
    }
}



static
void
get_exemplar_candidate_alignments(const starling_options& opt,
                                  const starling_deriv_options& dopt,
                                  const read_segment& rseg,
                                  const indel_synchronizer& isync,
                                  const alignment& exemplar,
                                  const known_pos_range realign_pr,
                                  mca_warnings& warn,
                                  std::set<candidate_alignment>& cal_set) {

    const unsigned read_length(rseg.read_size());

    starling_align_indel_status indel_status_map;
    std::vector<indel_key> indel_order;

    candidate_alignment cal;
    load_cal_with_edge_indels(exemplar,cal);

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
    // alignment. Note this won't work for grouper edge indels but we've
    // removed those from the exemplars already.
    //
    {
        indel_set_t cal_indels;
        get_alignment_indels(cal,opt.max_indel_size,cal_indels);

        BOOST_FOREACH(const indel_key& ik, cal_indels) {
            if(indel_status_map.find(ik)==indel_status_map.end()) {
                std::ostringstream oss;
                oss << "ERROR: Exemplar alignment contains indel not found in the overlap indel set\n"
                    << "\tIndel: " << ik
                    << "Exemplar overlap set:\n";
                dump_indel_status(indel_status_map,oss);
                throw blt_exception(oss.str().c_str());
            }
            indel_status_map[ik].is_present = true;
        }
    }

    // to prevent incompatible alignments, we must put all indels present in the exemplar first in the order list:
    //
    {
        typedef starling_align_indel_status::const_iterator siter;
        const siter i_begin(indel_status_map.begin());
        const siter i_end(indel_status_map.end());

        indel_order.clear();
        for(siter i(i_begin); i!=i_end; ++i) { if(i->second.is_present) indel_order.push_back(i->first); }
        for(siter i(i_begin); i!=i_end; ++i) { if(! i->second.is_present) indel_order.push_back(i->first); }
    }

    // to prevent truncated search, we must put all non-present remove-only indels last in the order list:
    sort_remove_only_indels_last(indel_status_map,indel_order);

#ifdef DEBUG_ALIGN
    std::cerr << "VARMIT exemplar starting indel_status_map:\n";
    dump_indel_status(indel_status_map,std::cerr);

    {
        std::cerr << "VARMIT exemplar starting indel_order:\n";
        const unsigned foo(indel_order.size());
        for(unsigned j(0); j<foo; ++j) {
            std::cerr << "no: " << j << " " << indel_order[j] << "\n";
        }
    }
#endif

    // to handle soft-clip and hard-clip in the genomic alignment, we take the
    // soft and hard clip portion off of the alignment before entering
    // make_candidate_alignments, and add it back in afterwards.
    //
    // TODO -- something less hacky to handle clipping
    //
    unsigned cal_read_length(read_length);
    unsigned hc_lead(0);
    unsigned hc_trail(0);
    unsigned sc_lead(0);
    unsigned sc_trail(0);
    const bool is_exemplar_clip(is_clipped(cal.al.path));
    if(is_exemplar_clip) {
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
        BOOST_FOREACH(candidate_alignment ical, cal_set2) {
            apath_clip_adder(ical.al.path,
                             hc_lead,
                             hc_trail,
                             sc_lead,
                             sc_trail);
            cal_set.insert(ical);
        }
    }

    {   // clear out-of-range alignment candidates:
        std::set<candidate_alignment> cal_set2(cal_set);
        cal_set.clear();
        BOOST_FOREACH(const candidate_alignment& ical, cal_set2) {
            // check that the alignment is within realign bounds
            if(is_alignment_spanned_by_range(realign_pr,ical.al)) {
                cal_set.insert(ical);
            }
        }
    }
}



// search for optimal realignments of the read and score alternate
// indel states in preparation for indel genotype calling
//
void
realign_and_score_read(const starling_options& opt,
                       const starling_deriv_options& dopt,
                       const starling_sample_options& sample_opt,
                       const reference_contig_segment& ref,
                       const known_pos_range& realign_pr,
                       read_segment& rseg,
                       indel_synchronizer& isync) {

    if(! rseg.is_valid()) {
        log_os << "ERROR: invalid alignment path associated with read segment:\n" << rseg;
        exit(EXIT_FAILURE);
    }

    // check that the original alignment is within realign bounds
    if(! is_alignment_spanned_by_range(realign_pr,rseg.genome_align())) return;

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
    get_exemplar_alignments(opt,rseg,exemplars);

    if(exemplars.empty()) return;

    // If any exemplars map to a negative start position, then skip realignment
    //
    // Note that even though such alignments will not occur in BAM, they
    // can still be produced by grouper.
    //
    BOOST_FOREACH(const alignment& al, exemplars) {
        if(al.pos<0) return;
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

    BOOST_FOREACH(const alignment& al, exemplars) {
        get_exemplar_candidate_alignments(opt,dopt,rseg,isync,al,realign_pr,warn,cal_set);
    }

    assert(! cal_set.empty());

    const bool is_incomplete_search(warn.origin_skip || warn.max_toggle_depth);

    // the max_toggle event is too common in genomic resequencing to have a
    // default warning:
    //
    const bool is_max_toggle_warn_enabled(opt.verbosity >= LOG_LEVEL::ALLWARN);

    if(warn.origin_skip || (warn.max_toggle_depth && is_max_toggle_warn_enabled)) {
        static const char* osr="alignments crossed chromosome origin";
        static const char* mir="exceeded max number of indel switches";
        const char* reason(warn.origin_skip ? osr : mir );


        log_os << "WARNING: re-alignment skipped some alternate alignments for read: "  << rseg.key()
               //               << " at buffer_pos: " << (sread.buffer_pos+1) << "\n"
               << "\treason: " << reason << "\n";
    }

    score_candidate_alignments_and_indels(opt,dopt,sample_opt,ref,rseg,isync,cal_set,is_incomplete_search);
}
