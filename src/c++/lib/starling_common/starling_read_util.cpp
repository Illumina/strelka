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

#include "alignment_util.hh"
#include "starling_read_util.hh"



pos_t
get_alignment_buffer_pos(const alignment& al)
{

    const pos_t lead(unalignedPrefixSize(al.path));
    return al.pos-lead;
}



namespace
{

struct ddata
{

    ddata(const unsigned read_size,
          const unsigned flank_size,
          read_mismatch_info& rmi_init)
        : is_totaled(false),
          fs(flank_size), fs2(fs*2),
          delta_size(std::max(1+fs2,read_size)-fs2),
          rmi(rmi_init)
    {
        for (unsigned i(0); i<delta_size; ++i)
        {
            rmi[i].delta = 0;
        }
    }

    void
    inc(const unsigned start_pos,
        const unsigned length)
    {
        assert(! is_totaled);
        rmi[std::max(fs2,start_pos)-fs2].delta += 1;
        if ((start_pos+length)<delta_size)
        {
            rmi[start_pos+length].delta -= 1;
        }
    }

    int
    get(const unsigned pos)
    {
        if (! is_totaled)
        {
            total();
            is_totaled=true;
        }
        return rmi[std::min(delta_size-1,std::max(fs,pos)-fs)].delta;
    }

private:
    void
    total()
    {
        for (unsigned i(1); i<delta_size; ++i)
        {
            rmi[i].delta += rmi[i-1].delta;
        }
    }

    bool is_totaled;
    const unsigned fs;
    const unsigned fs2;
    const unsigned delta_size;
    read_mismatch_info& rmi;
};

}



//
// basic algo (described with 1-based indices):
//
// S=length(seq) F=length(flank)
//
// 1. delta list length is DL=max(1,S-2F)
// 2. for each position:
//	1. if position is mismatch, add 1 to delta list at position max(1,p-2F), and subtract 1 from delta list at position p+1 if p+1 <= DL
// 3. sum delta list to match-count list
// 4. match(pos) = delta[min(DL,max(1,pos-F))]
//
//..for window size 1:
//
// match:       M
// pos  : 1  2  3  4  5
// delta: x  1  0  0  x
//
void
create_mismatch_filter_map(const blt_options& client_opt,
                           const alignment& al,
                           const bam_seq_base& ref_seq,
                           const bam_seq_base& read_seq,
                           const unsigned read_begin,
                           const unsigned read_end,
                           const CandidateSnvBuffer& candidateSnvBuffer,
                           read_mismatch_info& rmi)
{

    assert(read_begin<=read_end);
    const unsigned read_size(read_seq.size());
    assert(read_end<=read_size);
    const unsigned fs(client_opt.mismatchDensityFilterFlankSize);
    ddata dd(read_size,fs,rmi);

    for (unsigned i(0); i<read_size; ++i) rmi[i].is_mismatch=false;

    pos_t ref_head_pos(al.pos);
    unsigned read_head_pos(0);

    using namespace ALIGNPATH;

    const std::pair<unsigned,unsigned> ends(get_match_edge_segments(al.path));

    const unsigned as(al.path.size());
    for (unsigned i(0); i<as; ++i)
    {
        const path_segment& ps(al.path[i]);

        const bool is_edge_segment((i<ends.first) || (i>ends.second));

        if       (ps.type == INSERT)
        {
            if (! is_edge_segment) dd.inc(read_head_pos,ps.length);
            read_head_pos += ps.length;

        }
        else if (ps.type == DELETE)
        {
            if (! is_edge_segment) dd.inc(read_head_pos,0);
            ref_head_pos += ps.length;

        }
        else if (is_segment_align_match(ps.type))
        {
            for (unsigned j(0); j<ps.length; ++j)
            {
                const unsigned read_pos(read_head_pos+j);
                if ((read_pos < read_begin) || (read_pos >= read_end)) continue; // allow for read end trimming
                const pos_t ref_pos(ref_head_pos+static_cast<pos_t>(j));

                char readChar(read_seq.get_char(read_pos));
                if (readChar != ref_seq.get_char(ref_pos))
                {
                    // if the mismatch is a SNV found in an active region, don't increase the counter
                    if (not candidateSnvBuffer.isCandidateSnvAnySample(ref_pos, readChar))
                    {
                        rmi[read_pos].is_mismatch=true;
                        dd.inc(read_pos,1);
                    }
                }
            }
            read_head_pos += ps.length;
            ref_head_pos += ps.length;
        }
        else if (ps.type == SOFT_CLIP)
        {
            //dd.inc(read_head_pos,ps.length);
            read_head_pos += ps.length;

        }
        else if (ps.type == HARD_CLIP)
        {
            // do nothing
        }
        else
        {
            std::ostringstream oss;
            oss << "Can't handle cigar code: " << segment_type_to_cigar_code(ps.type) << "\n";
            throw blt_exception(oss.str().c_str());
        }
    }

    // set mismatch filter value:
    const int max_pass(static_cast<int>(client_opt.mismatchDensityFilterMaxMismatchCount));
    for (unsigned i(0); i<read_size; ++i)
    {
        const int del(dd.get(i));
        rmi[i].mismatch_count = del;
        rmi[i].mismatch_count_ns = del - rmi[i].is_mismatch;
        rmi[i].mismatch_filter_map = (max_pass < del);
    }
}



void
get_valid_alignment_range(const alignment& al,
                          const bam_seq_base& ref_seq,
                          const bam_seq_base& read_seq,
                          pos_range& valid_pr)
{
    static const int match_score(2);
    static const int mismatch_score(-5);
    static const int min_segment_score(-11);

    const unsigned read_size(read_seq.size());

    std::vector<int> fwd_read_score(read_size,0);
    std::vector<int> rev_read_score(read_size,0);

    pos_t ref_head_pos(al.pos);
    unsigned read_head_pos(0);

    using namespace ALIGNPATH;

    for (const path_segment& ps : al.path)
    {
        if       ((ps.type == INSERT) || (ps.type == SOFT_CLIP))
        {
            if (ps.type == INSERT)
            {
                fwd_read_score[read_head_pos] += mismatch_score;
                rev_read_score[read_head_pos+ps.length-1] += mismatch_score;
            }
            read_head_pos += ps.length;

        }
        else if (ps.type == DELETE)
        {

            if (read_head_pos > 0)   // filter out leading deletions:
            {
                fwd_read_score[read_head_pos-1] += mismatch_score;
            }

            if (read_head_pos < read_size)   // filter out trailing deletions
            {
                rev_read_score[read_head_pos] += mismatch_score;
            }
            ref_head_pos += ps.length;

        }
        else if (is_segment_align_match(ps.type))
        {
            for (unsigned j(0); j<ps.length; ++j)
            {
                const unsigned read_pos(read_head_pos+j);
                const pos_t ref_pos(ref_head_pos+static_cast<pos_t>(j));

                const char read_char(read_seq.get_char(read_pos));
                const char ref_char(ref_seq.get_char(ref_pos));
                if ((read_char != 'N') && (ref_char != 'N'))
                {
                    if (read_char != ref_char)
                    {
                        fwd_read_score[read_pos] += mismatch_score;
                        rev_read_score[read_pos] += mismatch_score;
                    }
                    else
                    {
                        fwd_read_score[read_pos] += match_score;
                        rev_read_score[read_pos] += match_score;
                    }
                }
            }
            read_head_pos += ps.length;
            ref_head_pos += ps.length;

        }
        else if (ps.type == HARD_CLIP)
        {
            // do nothing:
        }
        else
        {
            std::ostringstream oss;
            oss << "Can't handle cigar code: " << segment_type_to_cigar_code(ps.type) << "\n";
            throw blt_exception(oss.str().c_str());
        }
    }

    // The reckoning:
    valid_pr.set_begin_pos(0);
    valid_pr.set_end_pos(read_size);
    int fwd_sum(0),fwd_min(min_segment_score);
    int rev_sum(0),rev_min(min_segment_score);
    for (unsigned i(0); i<read_size; ++i)
    {
        fwd_sum += fwd_read_score[i];
        if (fwd_sum<=fwd_min)
        {
            valid_pr.begin_pos=i+1;
            fwd_min=fwd_sum;
        }
        rev_sum += rev_read_score[read_size-i-1];
        if (rev_sum<=rev_min)
        {
            valid_pr.end_pos=read_size-i-1;
            rev_min=rev_sum;
        }
    }
    if (valid_pr.end_pos<=valid_pr.begin_pos)
    {
        valid_pr.begin_pos=0;
        valid_pr.end_pos=0;
    }
}
