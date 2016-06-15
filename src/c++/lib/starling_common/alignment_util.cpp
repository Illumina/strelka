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

#include "alignment_util.hh"
#include "blt_util/align_path.hh"
#include "blt_util/align_path_util.hh"
#include "htsapi/align_path_bam_util.hh"

#include <cassert>

#include <algorithm>


known_pos_range
get_strict_alignment_range(const alignment& al)
{
    const pos_t asize(apath_ref_length(al.path));

    return known_pos_range(al.pos,al.pos+asize);
}



known_pos_range
get_soft_clip_alignment_range(const alignment& al)
{
    const pos_t lead(apath_insert_lead_size(al.path));
    const pos_t trail(apath_insert_trail_size(al.path));
    const pos_t asize(apath_ref_length(al.path));

    const pos_t begin_pos(al.pos-lead);
    pos_t end_pos(al.pos+asize+trail);

    return known_pos_range(begin_pos,end_pos);
}



known_pos_range
get_alignment_range(
    const alignment& al)
{
    const pos_t lead(apath_read_lead_size(al.path));
    const pos_t trail(apath_read_trail_size(al.path));
    const pos_t asize(apath_ref_length(al.path));

    const pos_t begin_pos(al.pos-lead);
    pos_t end_pos(al.pos+asize+trail);

    return known_pos_range(begin_pos,end_pos);
}



known_pos_range
get_alignment_zone(const alignment& al,
                   const unsigned seq_length)
{
    const known_pos_range ps(get_alignment_range(al));
    known_pos_range ps2(ps);
    ps2.begin_pos=std::max(0,std::min(ps.begin_pos,ps.end_pos-static_cast<pos_t>(seq_length)));
    ps2.end_pos=std::max(ps.end_pos,ps.begin_pos+static_cast<pos_t>(seq_length));

    return ps2;
}



bool
is_indel_in_alignment(const alignment& al,
                      const IndelKey& indelKey,
                      pos_range& read_indel_pr)
{
    using namespace ALIGNPATH;

    read_indel_pr.clear();

    const path_t& path(al.path);
    unsigned path_index(0);
    unsigned read_offset(0);
    pos_t ref_head_pos(al.pos);
    const std::pair<unsigned,unsigned> ends(get_match_edge_segments(al.path));
    const unsigned aps(path.size());
    while (path_index<aps)
    {
        if (ref_head_pos > indelKey.right_pos()) return false;

        const path_segment& ps(path[path_index]);

        const bool is_edge_segment((path_index<ends.first) || (path_index>ends.second));

        const bool is_swap_start(is_segment_swap_start(path,path_index));

        assert(! (ps.type == SKIP));
        assert(! (is_edge_segment && is_swap_start));

        unsigned n_seg(1); // number of path_segments consumed
        if       (is_edge_segment)
        {
            // only edge segment we use is insert:
            if (ps.type == INSERT)
            {
                if (path_index<ends.first)
                {
                    if (      (ref_head_pos==indelKey.pos) &&
                              (INDEL::BP_RIGHT==indelKey.type))
                    {
                        read_indel_pr.set_end_pos(read_offset+ps.length);
                        return true;
                    }
                    else if ((ref_head_pos==indelKey.right_pos()) &&
                             ((INDEL::INSERT==indelKey.type) || (INDEL::SWAP==indelKey.type)) &&
                             (ps.length <= indelKey.length))
                    {
                        read_indel_pr.set_end_pos(read_offset+ps.length);
                        return true;
                    }
                }

            }
            else
            {
                if (      (ref_head_pos==indelKey.pos) &&
                          (INDEL::BP_LEFT==indelKey.type))
                {
                    read_indel_pr.set_begin_pos(read_offset);
                    return true;
                }
                else if ((ref_head_pos==indelKey.pos) &&
                         ((INDEL::INSERT==indelKey.type) || (INDEL::SWAP==indelKey.type)) &&
                         (ps.length <= indelKey.length))
                {
                    read_indel_pr.set_begin_pos(read_offset);
                    return true;
                }
            }

        }
        else if (is_swap_start)
        {
            const swap_info sinfo(path,path_index);
            n_seg=sinfo.n_seg;

            if ((ref_head_pos==indelKey.pos) &&
                (sinfo.insert_length==indelKey.length) &&
                (sinfo.delete_length==indelKey.swap_dlength))
            {
                read_indel_pr.set_begin_pos(read_offset);
                read_indel_pr.set_end_pos(read_offset+sinfo.insert_length);
                return true;
            }

        }
        else if (is_segment_type_indel(path[path_index].type))
        {
            if ((ref_head_pos==indelKey.pos) &&
                (ps.length == indelKey.length) &&
                (((INSERT==ps.type) && (INDEL::INSERT==indelKey.type)) ||
                 ((DELETE==ps.type) && (INDEL::DELETE==indelKey.type))))
            {
                read_indel_pr.set_begin_pos(read_offset);
                const unsigned insert_length(INDEL::INSERT==indelKey.type ? ps.length : 0);
                read_indel_pr.set_end_pos(read_offset+insert_length);
                return true;
            }
        }

        for (unsigned i(0); i<n_seg; ++i)
        {
            increment_path(path,path_index,read_offset,ref_head_pos);
        }
    }

    return false;
}



alignment
remove_edge_deletions(
    const alignment& al,
    const bool is_remove_leading_edge,
    const bool is_remove_trailing_edge)
{
    using namespace ALIGNPATH;

    alignment al2;
    al2.is_fwd_strand=al.is_fwd_strand;
    al2.pos=al.pos;

    const std::pair<unsigned,unsigned> ends(get_match_edge_segments(al.path));
    const unsigned as(al.path.size());
    for (unsigned i(0); i<as; ++i)
    {
        const path_segment& ps(al.path[i]);
        const bool is_leading_edge_segment(i<ends.first);
        const bool is_trailing_edge_segment(i>ends.second);
        const bool is_target_type(ps.type==DELETE);
        if (is_target_type &&
            ((is_leading_edge_segment && is_remove_leading_edge) ||
             (is_trailing_edge_segment && is_remove_trailing_edge)))
        {
            if (i<ends.first) al2.pos += ps.length;
        }
        else
        {
            al2.path.push_back(ps);
        }
    }

    return al2;
}



/// convert segment_type to match if the segment exists before or after all match segments currently in the alignment
///
static
alignment
matchify_edge_segment_type(
    const alignment& al,
    const ALIGNPATH::align_t segment_type,
    const bool is_match_leading_edge = true,
    const bool is_match_trailing_edge = true)
{
    using namespace ALIGNPATH;

    assert(is_segment_type_read_length(segment_type));

    alignment al2;
    al2.is_fwd_strand=al.is_fwd_strand;
    al2.pos=al.pos;

    const std::pair<unsigned,unsigned> ends(get_match_edge_segments(al.path));
    const unsigned as(al.path.size());
    for (unsigned i(0); i<as; ++i)
    {
        const path_segment& ps(al.path[i]);
        const bool is_leading_edge_segment(i<ends.first);
        const bool is_trailing_edge_segment(i>ends.second);
        const bool is_target_type(ps.type==segment_type);
        const bool is_candidate_edge((is_match_leading_edge && is_leading_edge_segment) ||
                                     (is_match_trailing_edge && is_trailing_edge_segment));
        const bool is_edge_target(is_candidate_edge && is_target_type);
        if (is_edge_target && is_leading_edge_segment) al2.pos-=ps.length;
        if (is_edge_target || (is_segment_align_match(ps.type)))
        {
            if ((! al2.path.empty()) && (is_segment_align_match(al2.path.back().type)))
            {
                al2.path.back().length += ps.length;
            }
            else
            {
                al2.path.push_back(ps);
                al2.path.back().type = MATCH;
            }
        }
        else
        {
            al2.path.push_back(ps);
        }
    }

    return al2;
}



alignment
matchify_edge_insertions(
    const alignment& al,
    const bool is_match_leading_edge,
    const bool is_match_trailing_edge)
{
    return matchify_edge_segment_type(al, ALIGNPATH::INSERT,is_match_leading_edge,is_match_trailing_edge);
}



alignment
matchify_edge_indels(
    const alignment& al,
    const bool is_match_leading_edge,
    const bool is_match_trailing_edge)
{
    const alignment al2(remove_edge_deletions(al,is_match_leading_edge,is_match_trailing_edge));
    return matchify_edge_insertions(al2,is_match_leading_edge,is_match_trailing_edge);
}



alignment
matchify_edge_soft_clip(
    const alignment& al)
{
    return matchify_edge_segment_type(al, ALIGNPATH::SOFT_CLIP);
}



void
getAlignmentFromBamRecord(
    const bam_record& br,
    alignment& al)
{
    al.clear();
    al.pos=br.pos()-1;
    al.is_fwd_strand=br.is_fwd_strand();
    bam_cigar_to_apath(br.raw_cigar(),br.n_cigar(),al.path);
}



#if 0
bool
normalize_alignment(alignment& al,
                    const std::string& read_seq,
                    const std::string& ref_seq)
{
    assert(0);
}



// original version of code --very efficient for single indels but not designed to handle
// potential complexities of multiple indels interacting during normalization
//


// shift indels to occur as far "to the left" as possible
//
// return two pieces of info:
// 1) bool indicating if the indel is invalid (i.e. it 'fell off the left end of the read')
// 2) distance to shift the indel to the left
//
// Note that this functions purpose is not to validate indels -- that
// will happen in a subsequent step, if an indel 'falls off' the edge
// of a read during normalization, it's found to be invalid by
// happenstance, but there's no reason to check for the same evidence
// of a non-informative indel by trying to push it off the right side
// of the read
//
std::pair<bool,unsigned>
normalize_indel(const char* base_seq,
                const pos_t base_seq_start,
                const char* insert_seq,
                const pos_t insert_seq_start,
                const unsigned insert_size)
{

    assert(NULL != base_seq);
    assert(NULL != insert_seq);
    assert(0 != insert_size);

    unsigned bs(base_seq_start);
    unsigned is(insert_seq_start);

    while (true)
    {
        // read deletion fell off edge:
        if (bs==0) return std::make_pair(true,0);

        // attempt to shift back one base
        if (base_seq[bs-1] != insert_seq[is+insert_size-1])
        {
            return std::make_pair(false,base_seq_start-bs);
        }

        // read insertion fell off edge:
        if (is==0) return std::make_pair(true,0);
        bs--;
        is--;
    }
}



// this function is held up on normalization of multiple, potentially colliding
// indels -- really this is becoming just a bad way of doing a realignment. Punt
// this whole procedure for now.
//

// Shift all indels as far "to the left" as possible -- note that
// some indels may be lost. Returns true if the alignment was changed.
//
bool
normalize_alignment(alignment& al,
                    const std::string& read_seq,
                    const std::string& ref_seq)
{

    bool is_norm(false);

    std::pair<bool,unsigned> norm_res;
    unsigned read_offset(0);
    unsigned ref_offset(0);

    pos_t last_event_end(ps.pos);

    const unsigned aps(al.apath.size());
    for (unsigned i(0); i<aps; ++i)
    {
        path_segment& ps(al.apath[i]);

        std::cerr << "apath segment i: " << i << "\n";

        if ((ps.type == INSERT) or (ps.type == DELETE))
        {
            // do indel normalization:
            //
            const char* base(ref_seq.c_str());
            pos_t base_start(al.pos+ref_offset);
            const char* insert(al.seq());
            pos_t insert_start(read_offset);
            if (ps.type == DELETE)
            {
                std::swap(base,insert);
                std::swap(base_start,insert_start);
            }
            std::cerr << "norming\n";
            norm_res=normalize_indel(base,base_start,
                                     insert,insert_start,
                                     ps.length);

            if (norm_res.first or (0 != norm_res.second)) is_norm=true;

            if (norm.res_first)
            {
                if (ps.type == INSERT)
                {
                    ps.type == MATCH;
                    ps.pos -= ps.length;
                }
                else
                {
                    ps.type == NONE;
                    ps.pos += ps.length;
                }
            }
            else
            {
                const unsigned shift(norm_res.second); // TODO -- collision detection!
                std::cerr << "norm pass shift: " << shift << "\n";

                if (ps.length <= MAX_INDEL_SIZE)
                {
                    indel in;
                    in.pos=base_pos+ref_offset-shift;
                    in.key.length = ps.length;
                    if (ps.type == INSERT)
                    {
                        in.key.type=INDEL::INSERT;
                        in.data.seq=std::string(fwd_strand_read+read_offset-shift,ps.length);
                    }
                    else
                    {
                        in.key.type=INDEL::DELETE;
                    }
                    in.data.rid_set.insert(read_id);
                    sppr.insert_indel(in);
                }
                else
                {
                    // TODO setup large indel breakpoints:
                    log_os << "WARNING: skipping large indel\n";
                }
            }
        }

        if       (ps.type == MATCH)
        {
            read_offset += ps.length;
            ref_offset += ps.length;
            continue;
        }
        else if (ps.type == DELETE)
        {
            ref_offset += ps.length;
        }
        else if (ps.type == INSERT)
        {
            read_offset += ps.length;
        }
        else
        {
            assert(0); // can't handle other CIGAR types yet
        }

    }

    return is_norm;
}
#endif
