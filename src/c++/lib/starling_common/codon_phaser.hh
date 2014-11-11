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
/*
 * Codonphaser.hh
 *
 *  Test for codon-phasing.
 *
 *  Created on: Aug 10, 2013
 *  Author: Morten Kallberg
 */

#pragma once

#include "starling_common/gvcf_locus_info.hh"
#include "starling_common/starling_read_buffer.hh"
#include "starling_common/starling_shared.hh"
#include "starling_common/gvcf_locus_info.hh"

#include <climits>


struct Codon_phaser
{
    Codon_phaser(
        const starling_options& init_opt,
        const starling_deriv_options& init_dopt,
        starling_read_buffer& init_read_buffer,
        const unsigned init_max_read_len)
        : opt(init_opt),
          read_buffer(init_read_buffer),
          max_read_len(init_max_read_len)
    {
        this->last_cleared = init_dopt.report_range.begin_pos-this->max_read_len;
        block_start     = -1;
        phase_indels    = false;              //TODO not used; if false we break the block when encountering an indel
        last_cleared    = -1;
        this->clear_buffer();
    }

    /// add site to buffer
    ///
    /// \returns true when the buffer should be printed as a phased block
    bool add_site(const site_info& si);

    void clear_buffer();                // clear site buffer
    void write_out_buffer() const;      // debugging feature, print current buffer to std
    void write_out_alleles() const;     // print allele evidence

    /// Are we currently in a phasing block?
    bool is_in_block() const { return _is_in_block; }

    /// buffer of het snp calls
    ///
    const std::vector<site_info>&
    buffer() const
    {
        return _buffer;
    }

private:
    void make_record();                 // make phased record
    void clear_read_buffer(const int& pos);    // free up read that are no longer in phasing evidence, up to and including this position

    void
    collect_read_segment_evidence(
        const read_segment& rseg);

    void collect_read_evidence();       // fill in allele counter
    void construct_reference();         // assemble the reference allele for the record
    void create_phased_record();        // fill in the si record and decide if we have sufficient evidence for a phased call
    unsigned get_block_length() const
    {
        return (this->block_end-this->block_start+1);
    }

    bool _is_in_block;
    std::vector<site_info> _buffer;
    const starling_options& opt;
    starling_read_buffer& read_buffer;  // pass along the relevant read-buffer
    int max_read_len;                   // the length of the input reads

    int block_start,block_end;          // position of the first and last added het site to block
    int last_cleared;
    int het_count;                      // total hets observed in buffer
    int total_reads,total_reads_unused; // total used and unused reads spanning phasing region
    bool phase_indels;                  // should we attempt to phase indels as well? For now false, thus returning any block upon encountering an indel
    std::string reference;              // the phased allele reference
    typedef std::map<std::string,int> allele_map;
    allele_map observations;
};
