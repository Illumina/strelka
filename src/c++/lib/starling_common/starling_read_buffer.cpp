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
#include "CandidateAlignment.hh"
#include "starling_read_util.hh"

#include "starling_common/starling_read_buffer.hh"

#include <cassert>

#include <iostream>


const starling_read_buffer::segment_group_t starling_read_buffer::_empty_segment_group;



starling_read_buffer::
~starling_read_buffer()
{
    for (auto& val : _read_data) delete val.second;
}



align_id_t
starling_read_buffer::
add_read_alignment(
    const bam_record& br,
    const alignment& inputAlignment,
    const MAPLEVEL::index_t maplev)
{
    assert(! br.is_unmapped());

    const align_id_t readIndex(getNextReadIndex());
    _read_data[readIndex] = new starling_read(br, inputAlignment, maplev, readIndex);
    starling_read& sread(*(_read_data[readIndex]));

    if (sread.isSpliced())
    {
        // deal with spliced reads now:
        const uint8_t exonCount(sread.getExonCount());
        for (unsigned exonIndex(0); exonIndex<exonCount; ++exonIndex)
        {
            const uint8_t readSegmentIndex(exonIndex+1);
            auto& readSegment(sread.get_segment(readSegmentIndex));
            const pos_t seg_buffer_pos(get_alignment_buffer_pos(readSegment.getInputAlignment()));
            readSegment.buffer_pos = seg_buffer_pos;
            (_pos_group[seg_buffer_pos]).insert(std::make_pair(readIndex,readSegmentIndex));
        }
    }
    else
    {
        const pos_t buffer_pos(get_alignment_buffer_pos(inputAlignment));
        const seg_id_t seg_id(0);
        sread.get_full_segment().buffer_pos = buffer_pos;
        (_pos_group[buffer_pos]).insert(std::make_pair(readIndex,seg_id));
    }

    return readIndex;
}


#if 1
void
starling_read_buffer::
rebuffer_read_segment(const align_id_t read_id,
                      const seg_id_t seg_id,
                      const pos_t new_buffer_pos)
{
    // double check that the read exists:
    const read_data_t::iterator i(_read_data.find(read_id));
    if (i == _read_data.end()) return;

    read_segment& rseg(i->second->get_segment(seg_id));

    // remove from old pos list:
    const pos_group_t::iterator j(_pos_group.find(rseg.buffer_pos));
    assert(j != _pos_group.end());
    const segment_t segkey(std::make_pair(read_id,seg_id));
    assert(j->second.count(segkey)==1);
    j->second.erase(segkey);

    // alter data within read:
    rseg.buffer_pos=new_buffer_pos;

    // add to new pos list:
    (_pos_group[new_buffer_pos]).insert(segkey);
}
#endif



read_segment_iter
starling_read_buffer::
get_pos_read_segment_iter(const pos_t pos)
{
    const segment_group_t* g(&(_empty_segment_group));
    const pos_group_t::const_iterator j(_pos_group.find(pos));
    if (j != _pos_group.end()) g=(&(j->second));
    return read_segment_iter(*this,g->begin(),g->end());
}



void
starling_read_buffer::
clear_iter(
    const pos_group_t::iterator i)
{
    segment_group_t& seg_group(i->second);
    for (const auto& val : seg_group)
    {
        const align_id_t read_id(val.first);
        const seg_id_t seg_id(val.second);

        const read_data_t::iterator k(_read_data.find(read_id));
        if (k == _read_data.end()) continue;

        const starling_read* srp(k->second);

        // only remove read from data structure when we find the last
        // segment: -- note this assumes that two segments will not
        // occur at the same position:
        //
        if (seg_id != srp->getExonCount()) continue;

        // remove from simple lookup structures and delete read itself:
        _read_data.erase(k);

        delete srp;
    }
    _pos_group.erase(i);
}



void
starling_read_buffer::
dump_pos(const pos_t pos,
         std::ostream& os) const
{
    const pos_group_t::const_iterator i(_pos_group.find(pos));
    if (i == _pos_group.end()) return;

    os << "READ_BUFFER_POSITION: " << pos << " DUMP ON\n";

    const segment_group_t& seg_group(i->second);
    segment_group_t::const_iterator j(seg_group.begin()),j_end(seg_group.end());
    for (unsigned r(0); j!=j_end; ++j)
    {
        const align_id_t read_id(j->first);
        const seg_id_t seg_id(j->second);
        const read_data_t::const_iterator k(_read_data.find(read_id));
        if (k == _read_data.end()) continue;

        const starling_read& sr(*(k->second));
        os << "READ_BUFFER_POSITION: " << pos << " read_segment_no: " << ++r << " seg_id: " << seg_id << "\n";
        os << sr.get_segment(seg_id);
    }
    os << "READ_BUFFER_POSITION: " << pos << " DUMP OFF\n";
}



read_segment_iter::ret_val
read_segment_iter::
get_ptr()
{
    static const ret_val null_ret(std::make_pair(static_cast<starling_read*>(nullptr),0));
    if (_head==_end) return null_ret;
    const align_id_t read_id(_head->first);
    const seg_id_t seg_id(_head->second);
    const starling_read_buffer::read_data_t::iterator i(_buff._read_data.find(read_id));
    if (i==_buff._read_data.end()) return null_ret;
    return std::make_pair(i->second,seg_id);
}
