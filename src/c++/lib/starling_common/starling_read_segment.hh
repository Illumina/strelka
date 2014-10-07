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
    is_treated_as_tier1_mapping() const;

    // is this a tier1/tier2 mapping?:
    bool
    is_treated_as_anytier_mapping() const;

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
        return NULL;
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
