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
/// \author Chris Saunders
///

#include "starling_common/starling_read_segment.hh"
#include "starling_common/starling_read.hh"

#include <iostream>



bool
read_segment::
is_tier1_mapping() const
{
    return _sread.is_tier1_mapping();
}



bool
read_segment::
is_tier1or2_mapping() const
{
    return _sread.is_tier1or2_mapping();
}



unsigned
read_segment::
full_read_size() const
{
    return _sread._read_rec.read_size();
}


bam_seq
read_segment::
get_bam_read() const
{
    return bam_seq(bam_get_seq(_sread.get_brp()),_size,_offset);
}



const uint8_t*
read_segment::
qual() const
{
    return bam_get_qual(_sread.get_brp())+_offset;
}



align_id_t
read_segment::
getReadIndex() const
{
    return _sread.getReadIndex();
}



read_key
read_segment::
key() const
{
    return _sread.key();
}



MAPLEVEL::index_t
read_segment::
getInputAlignmentMapLevel() const
{
    return _sread._inputAlignmentMapLevel;
}



uint8_t
read_segment::
map_qual() const
{
    return _sread.map_qual();
}



std::pair<bool,bool>
read_segment::
get_segment_edge_pin() const
{
    std::pair<bool,bool> res(false,false);
    const seg_id_t n_seg(_sread.getExonCount());
    for (unsigned i(0); i<n_seg; ++i)
    {
        const seg_id_t seg_id(i+1);
        if (this==&(_sread.get_segment(seg_id)))
        {
            if (i!=0) res.first=true;
            if ((i+1)!=n_seg) res.second=true;
        }
    }
    return res;
}



// detect whether this read has any small (handlable) alignments:
bool
read_segment::
is_any_nonovermax(const unsigned max_indel_size) const
{
    return (! getInputAlignment().empty()) &&
           (! getInputAlignment().is_overmax(max_indel_size));
}



bool
read_segment::
is_valid() const
{
    const read_segment& rseg(*this);
    const alignment& al(rseg.getInputAlignment());

    if (al.empty()) return false;

    return (! (is_apath_invalid(al.path,rseg.read_size()) ||
               is_apath_starling_invalid(al.path)));
}



void
short_report(std::ostream& os,
             const read_segment& rseg)
{

    if (!rseg.getInputAlignment().empty()) os << "INPUT " << rseg.getInputAlignment();
    os << "is_realigned? " << rseg.is_realigned << "\n";
    if (rseg.is_realigned)
    {
        os << "REALIGN " << rseg.realignment;
    }
    os << "buffer_pos: " << rseg.buffer_pos << "\n";
}



// full report for read_segment is designed to be used independently
// of starling_read:
//
std::ostream&
operator<<(std::ostream& os,
           const read_segment& rseg)
{

    os << "key: " << rseg.key() << "\n";
    os << "id: " << rseg.getReadIndex() << "\n";

    const bam_seq bseq(rseg.get_bam_read());
    os << "seq:  " << bseq << "\n";
    os << "qual: ";
    {
        const unsigned rs(rseg.read_size());
        const uint8_t* qual(rseg.qual());
        for (unsigned i(0); i<rs; ++i) os << static_cast<char>(qual[i]+33);
    }
    os << "\n";

    short_report(os,rseg);

    return os;
}
