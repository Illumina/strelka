//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
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

#include "starling_common/starling_read_segment.hh"
#include "starling_common/starling_read.hh"

#include <iostream>



bool
read_segment::
is_tier1_mapping() const
{
    return sread().is_tier1_mapping();
}



bool
read_segment::
is_tier1or2_mapping() const
{
    return sread().is_tier1or2_mapping();
}



unsigned
read_segment::
full_read_size() const
{
    return sread()._read_rec.read_size();
}


bam_seq
read_segment::
get_bam_read() const
{
    return bam_seq(bam_get_seq(sread().get_brp()),_size,_offset);
}



const uint8_t*
read_segment::
qual() const
{
    return bam_get_qual(sread().get_brp())+_offset;
}



align_id_t
read_segment::
id() const
{
    return sread().id();
}



read_key
read_segment::
key() const
{
    return sread().key();
}



MAPLEVEL::index_t
read_segment::
genome_align_maplev() const
{
    return sread().genome_align_maplev;
}



uint8_t
read_segment::
map_qual() const
{
    return sread().map_qual();
}



bool
read_segment::
is_full_segment() const
{
    return (read_size() == sread()._read_rec.read_size());
}


std::pair<bool,bool>
read_segment::
get_segment_edge_pin() const
{
    std::pair<bool,bool> res(false,false);
    const seg_id_t n_seg(sread().segment_count());
    for (unsigned i(0); i<n_seg; ++i)
    {
        const seg_id_t seg_id(i+1);
        if (this==&(sread().get_segment(seg_id)))
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

    const read_segment& rseg(*this);

    if ((! rseg.genome_align().empty()) &&
        (! rseg.genome_align().is_overmax(max_indel_size))) return true;

    return false;
}



bool
read_segment::
is_valid() const
{
    const read_segment& rseg(*this);
    const alignment& al(rseg.genome_align());

    if (al.empty()) return false;

    return (! (is_apath_invalid(al.path,rseg.read_size()) ||
               is_apath_starling_invalid(al.path)));
}



void
short_report(std::ostream& os,
             const read_segment& rseg)
{

    if (! rseg.genome_align().empty()) os << "GENOME " << rseg.genome_align();
    os << "is_realigned? " << rseg.is_realigned << "\n";
    if (rseg.is_realigned)
    {
        //os << "REALIGN_path_log_lhood: " << rseg.realign_path_lnp << "\n";
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
    os << "id: " << rseg.id() << "\n";

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
