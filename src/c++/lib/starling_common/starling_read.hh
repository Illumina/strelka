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
#include "htsapi/bam_dumper.hh"
#include "starling_common/starling_read_key.hh"
#include "starling_common/starling_read_segment.hh"

#include "boost/utility.hpp"

#include <iosfwd>
#include <memory>


/// Represents a single read and associated per-read data as required by strelka calling models
///
/// This class mostly provides structure to handle rna-seq reads, in which case each exon is handled
/// semi-independently as a "read segment"
///
struct starling_read : private boost::noncopyable
{
    /// \param br a representation of the read's htslib BAM record
    /// \param inputAlignment read alignment proposed by a mapper or other external tool
    /// \param inputAlignmentMapLevel mapping confidence classification for the input mapping
    ///
    /// This ctor handles segment setup for spliced reads.
    ///
    /// Note that the alignment and map-level info can all be derived from the bam_record, but it is
    /// passed into the ctor here only becuase it would have had to been computed anyway given strelka's
    /// current worklow
    starling_read(
        const bam_record& br,
        const alignment& inputAlignment,
        const MAPLEVEL::index_t inputAlignmentMapLevel,
        const align_id_t readIndex);

    // This is not const because we update the BAM record with the best
    // alignment if the read has been realigned:
    void
    write_bam(bam_dumper& bamd);

    bool
    is_fwd_strand() const
    {
        return _read_rec.is_fwd_strand();
    }

    uint8_t map_qual() const
    {
        return _read_rec.map_qual();
    }

    MAPLEVEL::index_t
    getInputAlignmentMapLevel() const
    {
        return _inputAlignmentMapLevel;
    }

    align_id_t getReadIndex() const
    {
        return _readIndex;
    }

    read_key
    key() const
    {
        return read_key(_read_rec);
    }

    bool
    isSpliced() const
    {
        return (! _exonInfo.empty());
    }

    seg_id_t
    getExonCount() const
    {
        return _exonInfo.size();
    }

    /// \brief Request a segment of the read
    ///
    /// \param segmentIndex The index value of 0 is reserved for the full alignment, and subsequent values refer to
    ///                  each exon
    read_segment&
    get_segment(const seg_id_t segmentIndex)
    {
        if (segmentIndex > 0)
        {
            assert(isSpliced() && (segmentIndex <= getExonCount()));
            return _exonInfo[segmentIndex - 1];
        }
        return _full_read;
    }

    /// Const variant of get_segment method
    const read_segment&
    get_segment(const seg_id_t segmentIndex) const
    {
        if (segmentIndex > 0)
        {
            assert(isSpliced() && (segmentIndex <= getExonCount()));
            return _exonInfo[segmentIndex - 1];
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

    /// Update full read segment with a realignment which integrates all individual exon realignments
    void
    update_full_segment();

    /// Mapping quality category for the input read
    const MAPLEVEL::index_t _inputAlignmentMapLevel;

    /// Internal alignment index created and used only within Strelka
    const align_id_t _readIndex;
    bam_record _read_rec;
    read_segment _full_read;

    /// Store details of each exon if the read is spliced
    std::vector<read_segment> _exonInfo;
};


std::ostream& operator<<(std::ostream& os, const starling_read& sr);
