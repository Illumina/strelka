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
 * assembler.hh
 *
 *  Test for assembler.
 *
 *  Created on: Aug 10, 2014
 *  Author: Morten Kallberg
 */

#pragma once

#include "gvcf_locus_info.hh"
#include "site_info_stream.hh"
#include "starling_common/starling_base_shared.hh"
#include "starling_common/starling_read_buffer.hh"
#include "assembly/predictor.hh"
#include "applications/starling/gvcf_block_site_record.hh"
#include <climits>


struct assembly_streamer : public site_info_stream
{
	assembly_streamer(
        const starling_base_options& init_opt,
        const starling_deriv_options& init_dopt,
        starling_read_buffer& init_read_buffer,
        const unsigned init_max_read_len,
		const RegionTracker& nocompress_regions)
        : opt(init_opt),
          dopt(init_dopt),
          read_buffer(init_read_buffer),
          max_read_len(init_max_read_len),
          myPredictor(opt,dopt,nocompress_regions)
    {
        this->clear();
    }

    // implement site_info_stream methods
    bool add_site(site_info& si);
	bool add_indel(const pos_t pos,
				  const indel_key ik,
				  const starling_diploid_indel_core& dindel,
				  const starling_indel_report_info& iri,
				  const starling_indel_sample_report_info& isri){return true;}

	void flush(){} //TODO
    // clear all object data
    void clear();

    void write_out_buffer() const;      // debugging feature, print current buffer to std
    void write_out_alleles() const;     // print allele evidence

    /// Are we currently in a phasing block?
    bool is_in_block() const
    {
        return block_start != -1;
    }

    /// buffer of variant calls
    ///
    const std::deque<site_info>&
    buffer() const
    {
        return this->_site_buffer;
    }

    int block_start,block_end;                  // position of the first and last added variant site

private:
    void make_record();                 // make contig VCF record

    void
    collect_read_segment_evidence(
        const read_segment& rseg);

    void construct_reference();         // 1. Determine the reference allele for the record
    void collect_read_evidence();       // 2. Fill in allele counter based on assembled graph
    void assemble();                    // 3. Generate contigs from reads using some assembly method
    void rescore(std::stringstream &AD);// 4. Score based on reads realigned to contigs
    void create_contig_records();       // 5. Fill in the si record and decide if we have sufficient evidence for a phased call
    unsigned get_block_length() const
    {
        return (this->block_end-this->block_start+1);
    }

    bool keep_collecting();                 // keep the block going
    known_pos_range2 do_assemble();                     // Assemble the region that is currently in the buffer


    const starling_base_options& opt;
    const starling_deriv_options& dopt;
    starling_read_buffer& read_buffer;          // pass along the relevant read-buffer
    int max_read_len;                           // the length of the input reads

//    int block_start,block_end;                  // position of the first and last added variant site
    int var_count;                              // total variant count
    int total_reads,total_reads_unused;         // total used and unused reads spanning assembled region
    bool phase_indels = false;                  // should we attempt to phase indels as well? For now false, thus returning any block upon encountering an indel
    std::string reference;                      // the phased allele reference
    typedef std::map<std::string,int> allele_map;
    allele_map observations;
    known_pos_range2 asmRegion;
    predictor myPredictor;
    std::vector<std::string> asm_contigs;
};
