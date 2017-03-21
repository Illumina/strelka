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

/// \file
///
/// \author Chris Saunders
///

#pragma once

#include "blt_common/map_level.hh"
#include "starling_common/alignment.hh"
#include "starling_common/starling_read_key.hh"

#include <iosfwd>
#include <map>


typedef uint8_t seg_id_t;
struct starling_read;


// class spun-off from starling read to represent the notion of exons
// without forcing a re-rig the entire starling design: as of now
// read-segment captures the concept of a read (sub-)segment for the
// purposes of realignment/pileup, etc. while starling_read aggregates
// all information associated with the full cluster.
//
struct read_segment
{
    read_segment(const uint16_t size=0,
                 const uint16_t offset=0,
                 const starling_read* sread_ptr=nullptr)
        : is_realigned(false),
          is_invalid_realignment(false),
          buffer_pos(0),
          _size(size),
          _offset(offset),
          _sread_ptr(sread_ptr) {}

#if 0
    MAPLEVEL::index_t
    effective_maplevel() const;
#endif

    bool
    is_tier1_mapping() const;

    // is this a tier1/tier2 mapping?:
    bool
    is_tier1or2_mapping() const;

    unsigned read_size() const
    {
        return _size;
    }
    bam_seq get_bam_read() const;
    const uint8_t* qual() const;

    // are there any alignments without indels above max_indel_size?
    bool
    is_any_nonovermax(const unsigned max_indel_size) const;

    // check that information is valid and self-consistent;
    bool is_valid() const;

    /// original alignment before re-aligning
    const alignment&
    genome_align() const
    {
        return _genome_align;
    }

    // these methods just repeat from the parent read:
    align_id_t
    id() const;

    read_key
    key() const;

    MAPLEVEL::index_t
    genome_align_maplev() const;

    uint8_t
    map_qual() const;

    /// return a bool pair, indicating the bin status at the begin and end of the
    /// read segment, respectively
    ///
    /// A "pinned" end of the segment means that end cannot be realigned, typically
    /// because it is an exon bounded by gap segments to other exon sequences
    ///
    std::pair<bool,bool>
    get_segment_edge_pin() const;

    // returns NULL for non-realigned contig reads:
    const alignment*
    get_best_alignment() const
    {
        if       (is_realigned)
        {
            return &(realignment);
        }
        else if (! genome_align().empty())
        {
            return &(genome_align());
        }
        return nullptr;
    }

//    const bool is_mate_unmapped(){return this->sread().is_mate_unmapped();}

private:
    friend struct starling_read;

    bool
    is_full_segment() const;

    const starling_read& sread() const
    {
        return *_sread_ptr;
    }

private:
    alignment _genome_align;
public:
    alignment realignment;
    //    float realign_path_lnp;
    bool is_realigned;
    bool is_invalid_realignment;
    pos_t buffer_pos;

private:
    uint16_t _size;
    uint16_t _offset;
    const starling_read* _sread_ptr;
};


// debugging output:
void
short_report(std::ostream& os, const read_segment& rseg);

std::ostream& operator<<(std::ostream& os, const read_segment& rseg);
