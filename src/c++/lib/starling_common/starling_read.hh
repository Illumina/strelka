// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

/// \file
///
/// \author Chris Saunders
///

#pragma once

#include "blt_common/map_level.hh"
#include "blt_util/bam_dumper.hh"
#include "starling_common/starling_read_key.hh"
#include "starling_common/starling_read_segment.hh"
#include "starling_common/starling_shared.hh"

#include "boost/utility.hpp"

#include <iosfwd>
#include <memory>


namespace READ_ALIGN
{
enum index_t
{
    GENOME
};

const char*
label(const index_t i);
}



// helper class for starling_read to support the recently nailed-on
// notion of exons
//
struct starling_segmented_read
{
    starling_segmented_read(const seg_id_t size);

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
// captures the concept of read as required by starling: namely
// potentially aligned by both ELAND and/or GROUPER
//
// all alignments must be on the same strand and 'reasonably' proximate.
//
// all alignment info is fwd-strand
//
struct starling_read : private boost::noncopyable
{
    starling_read(const bam_record& br,
                  const bool is_bam_record_genomic);

    // bool
    // is_bam_record_genomic() {
    //     return _is_bam_record_genomic;
    // }

    void
    set_genomic_bam_record(const bam_record& br)
    {
        assert(! _is_bam_record_genomic);
        _read_rec = br;
        _is_bam_record_genomic=true;
    }

    // is new alignment compatible with pre-existing information?
    //
    bool
    is_compatible_alignment(const alignment& al,
                            const READ_ALIGN::index_t rat,
                            const starling_options& opt) const;

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


//    const bool
//    is_mate_unmapped(){return this->_read_rec.is_mate_unmapped();}

private:
    friend struct read_segment;

    const bam1_t*
    get_brp() const
    {
        return _read_rec.get_data();
    }

#if 0
    // returns as tier1 mapped if only a contig alignment exists
    MAPLEVEL::index_t
    effective_maplevel() const;
#endif

    bool
    is_treated_as_anytier_mapping() const;

    bool
    is_treated_as_tier1_mapping() const;

    // update full segment with sub-segment realignments
    void
    update_full_segment();


public:
    // mapping qualities of ELAND reads, does not apply to GROUPER:
    MAPLEVEL::index_t genome_align_maplev;

private:
    bool _is_bam_record_genomic; // indicates that this was the original (and thus, complete) alignment before grouper.
    align_id_t _id;
    bam_record _read_rec;
    read_segment _full_read;
    std::unique_ptr<starling_segmented_read> _segment_ptr;
};


std::ostream& operator<<(std::ostream& os, const starling_read& sr);
