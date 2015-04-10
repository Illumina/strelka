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
#include "assembly/contigAssembler.hh"
#include "assembly/dummy_assembler.hh"
#include "assembly/string_graph_assembler.hh"
#include "applications/starling/gvcf_block_site_record.hh"
#include <climits>


struct assembly_streamer : public site_info_stream
{
	assembly_streamer(
        const starling_base_options& init_opt,
        const starling_deriv_options& init_dopt,
        starling_read_buffer& init_read_buffer,
        const unsigned init_max_read_len,
		const RegionTracker& nocompress_regions
		)
        : opt(init_opt),
          dopt(init_dopt),
          read_buffer(init_read_buffer),
          max_read_len(init_max_read_len)
    {

		if (opt.assembly_model=="bed"){
			// log_os << "Initializing bed assembler" << std::endl;
			this->myAssembler = new dummy_assembler();
		}
		else if (opt.assembly_model=="string-graph"){
			//log_os << "Initializing default assembler" << std::endl;
			this->myAssembler = new string_graph_assembler();
		}
		else if (opt.assembly_model=="de-bruijn"){
			//log_os << "Initializing default assembler" << std::endl;
			this->myAssembler = new dummy_assembler();
		}
		else{
			//log_os << "Initializing default assembler" << std::endl;
			this->myAssembler = new dummy_assembler();
		}

		this->clear();
    }

    // implement site_info_stream methods
    // these buffer sites until flush() is called
    bool add_site(site_info& si);
    bool add_indel(const indel_info& ii);

    // assemble collected sites and output to consumer
    void flush();

    void write_out_buffer() const;      // debugging feature, print current buffer to std
    void write_out_alleles() const;     // print allele evidence

    // REMOVEME -- debug pass-through of contig info from predictor
    void add_contig(std::string const & contig)
    {
        asm_contigs.push_back(contig);
    }

protected:
    // clear all object data
    void clear();

    int block_start, block_end;         // position of the first and last added variant site

private:
    void make_record();                 // make contig VCF record

    void collect_read_segment_evidence(const read_segment& rseg);

    void construct_reference();         // 1. Determine the reference allele for the record
    void collect_read_evidence();       // 2. Fill in allele counter based on assembled graph
    void assemble();                    // 3. Generate contigs from reads using some assembly method
    void rescore(std::stringstream &AD);// 4. Score based on reads realigned to contigs
    void create_contig_records();       // 5. Fill in the si record and decide if we have sufficient evidence for a phased call

    const starling_base_options& opt;
    const starling_deriv_options& dopt;
    starling_read_buffer& read_buffer;          // pass along the relevant read-buffer
    int max_read_len;                           // the length of the input reads

    int var_count;                              // total variant count
    int total_reads,total_reads_unused;         // total used and unused reads spanning assembled region
    bool phase_indels = false;                  // should we attempt to phase indels as well? For now false, thus returning any block upon encountering an indel
    std::string reference;                      // the phased allele reference
    typedef std::map<std::string, int> allele_map;
    allele_map observations;
    known_pos_range2 asmRegion;
    predictor_interface myPredictor;
    contigAssembler* myAssembler;

    std::vector<std::string> asm_contigs;
};
