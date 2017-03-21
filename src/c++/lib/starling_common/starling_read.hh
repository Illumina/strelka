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

#pragma once

#include "blt_common/map_level.hh"
#include "htsapi/bam_dumper.hh"
#include "starling_common/starling_read_key.hh"
#include "starling_common/starling_read_segment.hh"

#include "boost/utility.hpp"

#include <iosfwd>
#include <memory>


// helper class for starling_read to support the recently nailed-on
// notion of exons
//
struct starling_segmented_read
{
    explicit
    starling_segmented_read(const seg_id_t size)
        : _seg_info(size) {}

    void
    set_segment(const seg_id_t seg_no,
                const read_segment& rseg);

    read_segment&
    get_segment(const seg_id_t seg_no);

    const read_segment&
    get_segment(const seg_id_t seg_no) const;

    seg_id_t
    segment_count() const
    {
        return _seg_info.size();
    }

private:
    std::vector<read_segment> _seg_info;
};



//
// captures the concept of read as required by starling
//
// all alignment info is fwd-strand
//
struct starling_read : private boost::noncopyable
{
    starling_read(const bam_record& br);

    // enters full alignment, and handles segment setup for splice
    // sites:
    void
    set_genome_align(const alignment& al);

    // nonconst because we update the BAM record with the best
    // alignment if the read has been realigned:
    void
    write_bam(bam_dumper& bamd);

    bool
    is_fwd_strand() const
    {
        return _read_rec.is_fwd_strand();
    }

    uint8_t map_qual() const;

    align_id_t& id()
    {
        return _id;
    }
    align_id_t id() const
    {
        return _id;
    }

    read_key
    key() const
    {
        return read_key(_read_rec);
    }

    bool
    is_segmented() const
    {
        return static_cast<bool>(_segment_ptr);
    }

    seg_id_t
    segment_count() const
    {
        if (is_segmented())
        {
            return _segment_ptr->segment_count();
        }
        return 0;
    }

    read_segment&
    get_segment(seg_id_t seg_no)
    {
        if (seg_no>0)
        {
            assert(is_segmented() && (seg_no<=segment_count()));
            return _segment_ptr->get_segment(seg_no);
        }
        return _full_read;
    }

    const read_segment&
    get_segment(seg_id_t seg_no) const
    {
        if (seg_no>0)
        {
            assert(is_segmented() && (seg_no<=segment_count()));
            return _segment_ptr->get_segment(seg_no);
        }
        return _full_read;
    }

    read_segment&
    get_full_segment()
    {
        return get_segment(0);
    }

    const read_segment&
    get_full_segment() const
    {
        return get_segment(0);
    }

private:
    friend struct read_segment;

    const bam1_t*
    get_brp() const
    {
        return _read_rec.get_data();
    }

    bool
    is_tier1_mapping() const;

    bool
    is_tier1or2_mapping() const;

    // update full segment with sub-segment realignments
    void
    update_full_segment();

public:
    // read mapper quality categories
    MAPLEVEL::index_t genome_align_maplev;

private:
    align_id_t _id;
    bam_record _read_rec;
    read_segment _full_read;
    std::unique_ptr<starling_segmented_read> _segment_ptr;
};


std::ostream& operator<<(std::ostream& os, const starling_read& sr);
