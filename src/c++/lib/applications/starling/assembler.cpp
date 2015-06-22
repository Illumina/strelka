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
 * assembler.cpp
 *
 *  Created on: Sep 10, 2014
 *  Author: Morten Kallberg
 */

#include <assembler.hh>
#include <array>
#include <sstream>
#include <vector>

#define DEBUG_ASSEMBLE

#ifdef DEBUG_ASSEMBLE
#include "blt_util/log.hh"
#endif

// Add a SNP site to the assembly buffer

bool
assembly_streamer::
add_site(site_info& si)
{
    this->_site_buffer.push_back(si);
    
    if(block_start < 0)
    {
        block_start = si.pos;
    }
    else
    {
        block_start = std::min(si.pos, block_start);
    }

    block_end = std::max(si.pos + 1, block_end);

    // Update block info
    if (si.dgt.is_snp)
    {
        var_count++;
    }

    return true;
}

bool
assembly_streamer::
add_indel(indel_info& ii)
{
    this->_indel_buffer.push_back(ii);

    if(block_start < 0)
    {
        block_start = ii.pos;
    }
    else
    {
        block_start = std::min(ii.pos, block_start);
    }

    //CASE: Start a new potential assembly block
    block_end = std::max(ii.pos + ii.ref_length() - 1, block_end);

    //std::cerr << "adding indel from " << ii.pos << " to " << block_end << " with length = " << ii.ref_length() << "\n";
    var_count++;

//    make_record();
//      log_os << "Assembling " << this->block_start << " - " << this->block_end << std::endl;
//    }
//    this->write_out_buffer();
//    this->notify_consumer();
//    this->clear();
    return true;
}

void
assembly_streamer::
flush()
{
    // modify gVCF records and buffer accordingly, notify consumer
    make_record();
//     log_os << "Assembling " << this->block_start << " - " << this->block_end << std::endl;
    clear();
}

// makes the phased VCF record from the buffered sites list
void
assembly_streamer::make_record()
{
    this->construct_reference();
//    this->collect_read_evidence();
    this->create_contig_records();
}

void
assembly_streamer::construct_reference()
{
    this->reference = "";
    // TODO: the following will need to be revised to handle indels!
    for (unsigned i=0; i<(_site_buffer.size()); i++)
    {
        if(_site_buffer.at(i).pos >= block_end)
        {
            break;
        }
        if(_site_buffer.at(i).pos >= block_start)
        {
            this->reference += _site_buffer.at(i).ref;
        }
    }
}

//void
//assembly_streamer::assemble()
//{
//    // THIS IS A STUB ASSEMBLER: IT ASSUMES THAT THE PREDICTOR CAN PROVIDE THE DESIRED ASSEMBLY
//    // The contigs have been passed in here already (hopefully) via the contig() function
//}

void
assembly_streamer::rescore(std::stringstream &AD)
{
    // TODO: some sort of read realignment-based scoring
    // TODO SOON: add score information here


    // add information from the alleles selected for output
    AD << 10;
//    bool hasAlt(false);
    for (const auto& val : asm_contigs)
    {
        if (val==reference) continue;
        AD << ',' << 15;
    }
    return;
}


void
assembly_streamer::create_contig_records()
{

    // pick first records in buffer as our anchoring point for the assembled record
    notify_consumer_up_to(block_start);
    site_info base = _site_buffer.at(0);
    site_info si = this->myAssembler->assembleReads(this->block_start,this->block_end,this->read_buffer,this->reference,base);
    clear_site_buffer_to_pos(block_end);
    clear_indel_buffer_to_pos(block_end);
    //assert("site buffer not empty" && _site_buffer.size() > 0);


    // add information from the alleles selected for output
//    std::stringstream alt;
//    bool hasAlt(false);
//    for (const auto& val : asm_contigs)
//    {
//        if (val==reference) continue;
//        hasAlt = true;
//        if (! alt.str().empty()) alt << ',';
//        alt << val;
//    }
//    if ( ! hasAlt )
//    {
//        alt << ".";
//    }

    std::stringstream AD;
    rescore(AD);

//    base.phased_ref = this->reference;
    //const bool is_ref(max_alleles[0].first==this->reference || max_alleles[1].first==this->reference);
    base.smod.is_block                  = false;
    base.smod.is_unknown                = false;
    // TODO: make this do reasonable things for GT
    base.smod.max_gt = 4;
    base.dgt.ref_gt  = 0; // setting the gt method to 0/1

    // set GQ and GQX - hard-coded for illustration
    base.smod.gq                = 30;
    base.dgt.genome.snp_qphred  = 1000;
    base.smod.gqx               = 30;
    base.Qscore                 = 20;
//    base.phased_alt = alt.str();
    base.phased_AD  = AD.str();

    // set more vcf record fields
//    int reads_ignored   = 10;
    base.n_used_calls   = 15;
    base.n_unused_calls = 20; // second term mark all alleles that we didnt use as unused reads

    //set any filters
//    base.smod.filters.reset();

    // specify that we are covering multiple records, and make the needed modification in the output
    base.smod.is_phased_region = true;
    base.smod.is_assembled_contig = true;

    // set more vcf record fields
//    int reads_ignored   = 10;
//    base.n_used_calls   = 15;
//    base.n_unused_calls = 20; // second term mark all alleles that we didnt use as unused reads

    // Add in assembled record(s)
    _site_buffer.push_front(base);

//    log_os << base << std::endl;
//    log_os << base.ref << std::endl;
    // write remaining site_infos
    this->notify_consumer_up_to(this->block_end);
    //    _site_buffer.pop_front();
}


void
assembly_streamer::
collect_read_evidence()
{
        // TODO hook in for assembler contigs
    this->observations[this->reference] = 10;
    std::string altAllele(this->reference.length(), 'N');
    this->observations[altAllele] = 15;
}

known_pos_range2
assembly_streamer::
do_assemble()
{
    //extend with more data structures to determine is assembly criteriea is met
    known_pos_range2 predAsmRegion; //= myPredictor.do_assemble_and_update();
    return predAsmRegion;
}


void
assembly_streamer::clear()
{
//    this->clear_buffer();
    observations.clear();
    block_start                 = -1;
    block_end                   = -1;
    var_count                   = 0;
    total_reads                 = 0;
    total_reads_unused          = 0;
    reference                   = "";
    asmRegion = known_pos_range2(-1,-1);
}

void
assembly_streamer::write_out_buffer() const
{
    for (const auto& val : _site_buffer)
    {
        log_os << val << " ref " << val.ref << "\n";
    }
}

void
assembly_streamer::write_out_alleles() const
{
    for (const auto& val : observations)
    {
        if (val.first == this->reference) continue;
        log_os << val.first << "=" << val.second << " ";
    }
    log_os << "\n";
}
