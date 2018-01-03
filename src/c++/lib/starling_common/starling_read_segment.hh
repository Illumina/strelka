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

#pragma once

#include "blt_common/map_level.hh"
#include "starling_common/alignment.hh"
#include "starling_common/starling_read_key.hh"

#include <iosfwd>


typedef uint8_t seg_id_t;
struct starling_read;


/// spin-off from starling read to generalize full-read alignments with exon alignments
///
/// read-segment captures the concept of a read (sub-)segment for the
/// purposes of realignment/pileup, etc. while starling_read aggregates
/// all information associated with the full read
///
struct read_segment
{
    /// \param size Length of the read segment
    /// \param offset Offset of the read segment within the full read
    /// \param sread Reference to the full read
    /// \param inputAlignment Alignment of the read segment before realignment (must not be empty)
    read_segment(
        const uint16_t size,
        const uint16_t offset,
        const starling_read& sread,
        const alignment& inputAlignment)
        : is_realigned(false),
          is_invalid_realignment(false),
          buffer_pos(0),
          _size(size),
          _offset(offset),
          _sread(sread),
          _inputAlignment(inputAlignment)
    {
        assert(! _inputAlignment.empty());
    }

    bool
    is_tier1_mapping() const;

    /// is this a tier1/tier2 mapping?:
    bool
    is_tier1or2_mapping() const;

    /// provides the size of the segment, not the full read
    unsigned read_size() const
    {
        return _size;
    }

    /// provides the size of the full read from which the segment is aligned
    ///
    /// this should always equal read_size except for RNA-Seq reads with GAPed alignments
    unsigned
    full_read_size() const;

    /// provides the offset of the segment in the full read
    ///
    /// this should always be zero except for RNA-Seq reads with GAPed alignments
    unsigned
    full_read_offset() const
    {
        return _offset;
    }

    bam_seq get_bam_read() const;
    const uint8_t* qual() const;

    /// are there any alignments without indels above max_indel_size?
    bool
    is_any_nonovermax(const unsigned max_indel_size) const;

    /// check that object is valid and self-consistent;
    bool is_valid() const;

    /// Get the read's input alignment prior to any realignment or clipping modifications
    const alignment&
    getInputAlignment() const
    {
        return _inputAlignment;
    }

    // these methods just repeat from the parent read:
    align_id_t
    getReadIndex() const;

    read_key
    key() const;

    MAPLEVEL::index_t
    getInputAlignmentMapLevel() const;

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

    /// Return either the read's starting alignment or its realignment if it has been realigned
    const alignment&
    getBestAlignment() const
    {
        if (is_realigned) return realignment;
        return getInputAlignment();
    }

    alignment realignment;
    bool is_realigned;
    bool is_invalid_realignment;
    pos_t buffer_pos;

private:
    uint16_t _size;
    uint16_t _offset;
    const starling_read& _sread;
    /// Read alignment as provided from the alignment input (BAM file, etc...)
    const alignment _inputAlignment;
};


// debugging output:
void
short_report(std::ostream& os, const read_segment& rseg);

std::ostream& operator<<(std::ostream& os, const read_segment& rseg);
