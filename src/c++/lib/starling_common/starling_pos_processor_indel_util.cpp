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
///
/// \author Chris Saunders
///
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the begining of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///



#include "starling_pos_processor_indel_util.hh"
#include "starling_read_util.hh"

#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"
#include "blt_util/align_path_util.hh"

#include <cassert>

#include <iostream>
#include <sstream>



static
void
finish_indel_sppr(IndelObservation& obs,
                  starling_pos_processor_base& posProcessor,
                  const unsigned sample_no)
{
    posProcessor.insert_indel(obs,sample_no);
}



static
void
bam_seq_to_str(
    const bam_seq_base& bs,
    const unsigned start,
    const unsigned end,
    std::string& s)
{
    s.clear();
    for (unsigned i(start); i<end; ++i) s.push_back(bs.get_char(i));
}



const unsigned max_cand_filter_insert_size(10);



/// Handle edge inserts. For contigs we take the edge insert verbatim. For contig
// reads we need to find out what (usually larger) indel this edge corresponds to.
/// for regular reads we record this as a private indel (ie. it does not count towards
/// candidacy.
///
static
void
process_edge_insert(
    const unsigned max_indel_size,
    const ALIGNPATH::path_t& path,
    const bam_seq_base& bseq,
    starling_pos_processor_base& posProcessor,
    IndelObservation& obs,
    const unsigned sample_no,
    const unsigned path_index,
    const unsigned read_offset,
    const pos_t ref_head_pos,
    const bool is_pinned_indel)
{
    using namespace ALIGNPATH;

    const path_segment& ps(path[path_index]);

    {
        // do not allow edge-indels on genomic reads to generate or support candidate
        // indels, except for pinned cases:
        if (! is_pinned_indel) return;
        if (ps.length > max_indel_size) return;

        obs.key.pos=ref_head_pos;
        assert(ps.type == INSERT);
        obs.key.type=INDEL::INDEL;

#ifdef SPI_DEBUG
        log_os << "FOOBAR: adding pinned edge indel: " << obs.key << "\n";
#endif

        bam_seq_to_str(bseq,read_offset,read_offset+ps.length,obs.key.insertSequence);

        finish_indel_sppr(obs,posProcessor,sample_no);
    }
}



/// handle edge deletions
///
/// currently these are only allowed for indels adjacent to an intron
///
static
void
process_edge_delete(
    const unsigned max_indel_size,
    const ALIGNPATH::path_t& path,
    starling_pos_processor_base& posProcessor,
    IndelObservation& obs,
    const unsigned sample_no,
    const unsigned path_index,
    const pos_t ref_head_pos,
    const bool is_pinned_indel)
{
    using namespace ALIGNPATH;

    const path_segment& ps(path[path_index]);

    // do not allow edge-indels on genomic reads to generate or support candidate
    // indels, except for pinned cases:
    if (! is_pinned_indel) return;
    if (ps.length > max_indel_size) return;

    obs.key.pos=ref_head_pos;
    obs.key.deletionLength = ps.length;
    assert(ps.type == DELETE);
    obs.key.type=INDEL::INDEL;

#ifdef SPI_DEBUG
    log_os << "FOOBAR: adding pinned edge indel: " << obs.key << "\n";
#endif
    finish_indel_sppr(obs,posProcessor,sample_no);
}



// Note that unlike other indel processors, the swap processor returns
// the number of segments consumed.
//
// Like regular indels, a swap will be reported as a breakpoint if
// it's too long (meaning longeset of insert,delete size).
//
static
unsigned
process_swap(
    const unsigned max_indel_size,
    const ALIGNPATH::path_t& path,
    const bam_seq_base& bseq,
    starling_pos_processor_base& posProcessor,
    IndelObservation& obs,
    const unsigned sample_no,
    const unsigned path_index,
    const unsigned read_offset,
    const pos_t ref_head_pos)
{
    using namespace ALIGNPATH;

    const swap_info sinfo(path,path_index);
    const unsigned swap_size(std::max(sinfo.insert_length,sinfo.delete_length));

    // large insertions are not filtered as noise:
    if (obs.data.is_noise)
    {
        if (sinfo.insert_length > max_cand_filter_insert_size)
        {
            obs.data.is_noise=false;
        }
    }

    if (swap_size <= max_indel_size)
    {
        obs.key.pos=ref_head_pos;
        obs.key.deletionLength=sinfo.delete_length;
        obs.key.type = INDEL::INDEL;
        bam_seq_to_str(bseq,read_offset,read_offset+sinfo.insert_length,obs.key.insertSequence);
        finish_indel_sppr(obs,posProcessor,sample_no);

    }
    else
    {
        // left side BP:
        {
            obs.key.pos=ref_head_pos;
            obs.key.type=INDEL::BP_LEFT;
            const unsigned start(read_offset);
            const unsigned size(bseq.size()-read_offset);
            const unsigned end(start+std::min(size,max_indel_size));
            bam_seq_to_str(bseq,start,end,obs.data.breakpointInsertionSequence);
            finish_indel_sppr(obs,posProcessor,sample_no);
        }

        // right side BP:
        {
            obs.key.pos=ref_head_pos+sinfo.delete_length;
            obs.key.type=INDEL::BP_RIGHT;
            const unsigned next_read_offset(read_offset+sinfo.insert_length);
            const unsigned start_offset(next_read_offset-std::min(next_read_offset,max_indel_size));
            bam_seq_to_str(bseq,start_offset,next_read_offset,obs.data.breakpointInsertionSequence);
            finish_indel_sppr(obs,posProcessor,sample_no);
        }
    }

    return sinfo.n_seg;
}



/// Handle the regular ol' insertions and deletions. Reports these
/// types as breakpoints when they're too long:
///
static
void
process_simple_indel(
    const unsigned max_indel_size,
    const ALIGNPATH::path_t& path,
    const bam_seq_base& bseq,
    starling_pos_processor_base& posProcessor,
    IndelObservation& obs,
    const unsigned sample_no,
    const unsigned path_index,
    const unsigned read_offset,
    const pos_t ref_head_pos)
{
    using namespace ALIGNPATH;

    const path_segment& ps(path[path_index]);

    // large insertion breakpoints are not filtered as noise:
    if (obs.data.is_noise)
    {
        if ((ps.type == INSERT) &&
            (ps.length > max_cand_filter_insert_size))
        {
            obs.data.is_noise=false;
        }
    }

    if (ps.length <= max_indel_size)
    {
        obs.key.pos=ref_head_pos;
        obs.key.type=INDEL::INDEL;
        if (ps.type == DELETE)
        {
            obs.key.deletionLength = ps.length;
        }
        else
        {
            bam_seq_to_str(bseq,read_offset,read_offset+ps.length,obs.key.insertSequence);
        }
        finish_indel_sppr(obs,posProcessor,sample_no);
    }
    else
    {
        // left side BP:
        {
            obs.key.pos=ref_head_pos;
            obs.key.type=INDEL::BP_LEFT;
            const unsigned start(read_offset);
            const unsigned size(bseq.size()-read_offset);
            const unsigned end(start+std::min(size,max_indel_size));
            bam_seq_to_str(bseq,start,end,obs.data.breakpointInsertionSequence);
            finish_indel_sppr(obs,posProcessor,sample_no);
        }
        // right side BP:
        {
            obs.key.pos=ref_head_pos;
            if (ps.type == DELETE) obs.key.pos+=ps.length;
            obs.key.type=INDEL::BP_RIGHT;

            const unsigned next_read_offset(read_offset+((ps.type==INSERT) ? ps.length : 0));
            const unsigned start_offset(next_read_offset-std::min(next_read_offset,max_indel_size));
            bam_seq_to_str(bseq,start_offset,next_read_offset,obs.data.breakpointInsertionSequence);
            finish_indel_sppr(obs,posProcessor,sample_no);
        }
    }
}



unsigned
addAlignmentIndelsToPosProcessor(
    const unsigned max_indel_size,
    const reference_contig_segment& ref,
    const alignment& al,
    const bam_seq_base& read_seq,
    starling_pos_processor_base& posProcessor,
    const INDEL_ALIGN_TYPE::index_t iat,
    const align_id_t id,
    const unsigned sample_no,
    const std::pair<bool, bool>& edge_pin,
    const bool isLowMapQuality)
{
    using namespace ALIGNPATH;

    const unsigned seq_len(read_seq.size());
    if (is_apath_invalid(al.path,seq_len))
    {
        std::ostringstream oss;
        oss << "Can't handle alignment path '" << apath_to_cigar(al.path) << "' -- " << get_apath_invalid_reason(al.path,seq_len);
        throw blt_exception(oss.str().c_str());
    }

    if (is_apath_starling_invalid(al.path))
    {
        std::ostringstream oss;
        oss << "Can't handle alignment path '" << apath_to_cigar(al.path) << "'";
        throw blt_exception(oss.str().c_str());
    }

    const rc_segment_bam_seq ref_bseq(ref);

    const std::pair<unsigned,unsigned> ends(get_match_edge_segments(al.path));

    pos_range valid_pr;
    get_valid_alignment_range(al,ref_bseq,read_seq,valid_pr);

    unsigned path_index(0);
    unsigned read_offset(0);
    pos_t ref_head_pos(al.pos);

    unsigned total_indel_ref_span_per_read(0);
    const unsigned aps(al.path.size());

    const bool isUsingActiveRegionDetector(posProcessor.is_active_region_detector_enabled());
    ActiveRegionReadBuffer& activeRegionReadBuffer(posProcessor.getActiveRegionReadBuffer(sample_no));
    if (isUsingActiveRegionDetector)
    {
        activeRegionReadBuffer.setAlignInfo(id, sample_no, iat, al.is_fwd_strand);
    }

    while (path_index<aps)
    {
        const path_segment& ps(al.path[path_index]);
        const bool is_begin_edge(path_index<ends.first);
        const bool is_end_edge(path_index>ends.second);
        const bool is_edge_segment(is_begin_edge || is_end_edge);

        const bool is_swap_start(is_segment_swap_start(al.path,path_index));

        assert(ps.type != SKIP);

        // fix STARKA-43:
        // assert(! (is_edge_segment && is_swap_start));

        IndelObservation obs;
        obs.data.iat = iat;
        obs.data.id = id;
        obs.data.is_low_map_quality = isLowMapQuality;

        if (! is_segment_align_match(ps.type))
        {
//            log_os << al.path <<" \n";
            pos_range indel_read_pr;
            indel_read_pr.set_begin_pos((read_offset==0) ? 0 : (read_offset-1));

            unsigned rlen(0);
            if       (is_swap_start)
            {
                const swap_info sinfo(al.path,path_index);
                rlen=sinfo.insert_length;

                if (sinfo.delete_length<=max_indel_size)
                {
                    total_indel_ref_span_per_read += sinfo.delete_length;
                }
            }
            else if (is_segment_type_read_length(ps.type))
            {
                rlen=ps.length;
            }
            else
            {
                if (ps.type == DELETE)
                {
                    if (ps.length <= max_indel_size)
                    {
                        total_indel_ref_span_per_read += ps.length;
                    }
                }
            }
            indel_read_pr.set_end_pos(std::min(seq_len,read_offset+1+rlen));
            if (! valid_pr.is_superset_of(indel_read_pr)) obs.data.is_noise=true;
        }

        unsigned n_seg(1); // number of path segments consumed
        if (is_edge_segment)
        {
            // is this indel occurring on a pinned edge (ie against an exon?)
            const bool is_pinned_indel((is_begin_edge && edge_pin.first) ||
                                       (is_end_edge && edge_pin.second));

            // edge inserts are allowed for intron adjacent and grouper reads, edge deletions for intron adjacent only:
            if (ps.type == INSERT)
            {
                process_edge_insert(max_indel_size,al.path,read_seq,
                                    posProcessor,obs,sample_no,
                                    path_index,read_offset,ref_head_pos,
                                    is_pinned_indel);
            }
            else if (ps.type == DELETE)
            {
                if (is_pinned_indel)
                {
                    process_edge_delete(max_indel_size, al.path, posProcessor, obs, sample_no, path_index, ref_head_pos,
                                        is_pinned_indel);
                }
            }
            else if (ps.type == SOFT_CLIP)
            {
                if (isUsingActiveRegionDetector && !isLowMapQuality)
                {
                    pos_t softClipStartPos(ref_head_pos);
                    // for leading/trailing soft-clip, place the segment at the end/start position
                    if (is_begin_edge)
                        --softClipStartPos;

                    const pos_t refPos(softClipStartPos);
                    std::string segmentSeq;
                    bam_seq_to_str(read_seq, read_offset, read_offset+ps.length, segmentSeq);
                    activeRegionReadBuffer.insertSoftClipSegment(id, refPos, segmentSeq, is_begin_edge);
                }
            }
        }
        else if (is_swap_start)
        {
            n_seg = process_swap(max_indel_size,al.path,read_seq,
                                 posProcessor,obs,sample_no,
                                 path_index,read_offset,ref_head_pos);

        }
        else if (is_segment_type_indel(al.path[path_index].type))
        {
//            log_os << read_offset << "\n";
//            log_os << int(rs.mapq) << "\n";
//            log_os << al.path << "\n";
//            log_os << int(rs.qual[read_offset]) << "\n\n";

            //obs.key.addRanksumInfo(static_cast<int>(rs.mapq),static_cast<int>(rs.qual[read_offset]),true); //TODO for updateing indel ranksums
            process_simple_indel(max_indel_size,al.path,read_seq,
                                 posProcessor,obs,sample_no,
                                 path_index,read_offset,ref_head_pos);
        }
        else if (isUsingActiveRegionDetector && !isLowMapQuality && is_segment_align_match(ps.type))
        {
            // detect active regions (match/mismatch)
            for (unsigned j(0); j < ps.length; ++j)
            {
                const pos_t ref_pos(ref_head_pos + static_cast<pos_t>(j));
                pos_t read_pos = read_offset + j;

                char base_char = read_seq.get_char(read_pos);

                if (ref.get_base(ref_pos) != base_char)
                {
                    activeRegionReadBuffer.insertMismatch(id, ref_pos, base_char);
                }
                else
                {
                    activeRegionReadBuffer.insertMatch(id, ref_pos);
                }
            }
        }

        for (unsigned i(0); i<n_seg; ++i)
        {
            increment_path(al.path,path_index,read_offset,ref_head_pos);
        }
    }

    return total_indel_ref_span_per_read;
}
