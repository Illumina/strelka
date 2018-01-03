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

#include "CandidateAlignment.hh"

#include "blt_util/align_path_util.hh"

#include <cassert>

#include <iostream>



std::ostream&
operator<<(
    std::ostream& os,
    const CandidateAlignment& cal)
{
    os << "CANDIDATE_ALIGNMENT: " << cal.al;
    os << "CANDIDATE_ALIGNMENT_LEAD: " << cal.leading_indel_key;
    os << "CANDIDATE_ALIGNMENT_TRAIL: " << cal.trailing_indel_key;

    return os;
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



void
getAlignmentIndels(
    const CandidateAlignment& cal,
    const reference_contig_segment& ref,
    const read_segment& rseg,
    const unsigned max_indel_size,
    const bool includeMismatches,
    indel_set_t& indels)
{
    using namespace ALIGNPATH;

    indels.clear();

    const path_t& path(cal.al.path);
    unsigned path_index(0);
    unsigned read_offset(0);
    pos_t ref_head_pos(cal.al.pos);

    const rc_segment_bam_seq refBamSeq(ref);

    const std::pair<unsigned,unsigned> ends(get_match_edge_segments(cal.al.path));
    const unsigned aps(path.size());
    while (path_index<aps)
    {
        const bool is_edge_segment((path_index<ends.first) || (path_index>ends.second));
        const bool is_swap_start(is_segment_swap_start(path,path_index));

        const path_segment& ps(path[path_index]);

        assert(! (ps.type == SKIP));

        unsigned n_seg(1); // number of path_segments consumed
        if (is_swap_start)
        {
            const swap_info sinfo(path,path_index);
            n_seg=sinfo.n_seg;
        }

        if       (is_edge_segment)
        {
            // ignore all edge segments except INSERT/DELETE (this includes SWAP):
            if ((DELETE == ps.type) ||
                (INSERT == ps.type))
            {
                if (path_index<ends.first)
                {
                    assert(cal.leading_indel_key.type != INDEL::NONE);
                    indels.insert(cal.leading_indel_key);
                }
                else
                {
                    assert(cal.trailing_indel_key.type != INDEL::NONE);
                    indels.insert(cal.trailing_indel_key);
                }
            }

        }
        else if (is_swap_start)
        {
            const swap_info sinfo(path,path_index);

            const unsigned swap_size(std::max(sinfo.insert_length,sinfo.delete_length));
            if (swap_size <= max_indel_size)
            {
                std::string insertSequence;
                bam_seq_to_str(rseg.get_bam_read(),read_offset,read_offset+sinfo.insert_length, insertSequence);
                indels.insert(IndelKey(ref_head_pos,INDEL::INDEL,sinfo.delete_length, insertSequence.c_str()));
            }
            else
            {
                indels.insert(IndelKey(ref_head_pos,INDEL::BP_LEFT));
                const pos_t right_pos(ref_head_pos+sinfo.delete_length);
                indels.insert(IndelKey(right_pos,INDEL::BP_RIGHT));
            }

        }
        else if (is_segment_type_indel(path[path_index].type))
        {
            if (ps.length <= max_indel_size)
            {
                IndelKey indelKey(ref_head_pos,INDEL::INDEL);
                if (INSERT == ps.type)
                {
                    bam_seq_to_str(rseg.get_bam_read(),read_offset,read_offset+ps.length,indelKey.insertSequence);
                }
                else
                {
                    indelKey.deletionLength = ps.length;
                }
                indels.insert(indelKey);
            }
            else
            {
                indels.insert(IndelKey(ref_head_pos,INDEL::BP_LEFT));
                const pos_t right_pos(ref_head_pos+((ps.type==INSERT) ? 0 : ps.length));
                indels.insert(IndelKey(right_pos,INDEL::BP_RIGHT));
            }

        }
        else if (includeMismatches && (is_segment_align_match(path[path_index].type)))
        {
            // match or mismatch
            for (unsigned i(0); i<ps.length; ++i)
            {
                const pos_t readPos(static_cast<pos_t>(read_offset + i));
                const uint8_t sbase(rseg.get_bam_read().get_code(readPos));
                if ((sbase == BAM_BASE::REF) || (sbase == BAM_BASE::ANY)) continue;
                const pos_t refPos(static_cast<pos_t>(ref_head_pos + i));
                bool isRef=(sbase == refBamSeq.get_code(refPos));
                if (isRef) continue;

                indels.insert(IndelKey(refPos, INDEL::MISMATCH, 1, std::string(1,get_bam_seq_char(sbase)).c_str()));
            }
        }
        for (unsigned i(0); i<n_seg; ++i)
        {
            increment_path(path,path_index,read_offset,ref_head_pos);
        }
    }
}
