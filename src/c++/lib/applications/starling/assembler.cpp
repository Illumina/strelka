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
 * Codon_phaser.cpp
 *
 *  Created on: Sep 10, 2013
 *  Author: Morten Kallberg
 */

#include <assembler.hh>
#include <array>
#include <sstream>
#include <vector>

//#define DEBUG_ASSEMBLE


#ifdef DEBUG_ASSEMBLE
#include "blt_util/log.hh"
#endif


// Add a SNP site to the phasing buffer
//,const gvcf_block_site_record& empty_block)

//bool add_indel(const pos_t pos,
//                        const indel_key ik,
//                        const starling_diploid_indel_core& dindel,
//                        const starling_indel_report_info& iri,
//                        const starling_indel_sample_report_info& isri){
//      return true;
//}

bool
assembly_streamer::
add_site(site_info& si)
{
//      return this->_consumer->add_site(si);

        this->_site_buffer.push_back(si);

        //CASE: Start a new potential assembly block
        block_end = si.pos+1;
    if (!is_in_block()){
        block_start = si.pos;
//        log_os << "New block case - starting @ " << si.pos  << std::endl;
    }

    // Update block info
    if (si.is_nonref())
        var_count ++;

    // CASE: Check if we should be extending the block further based on criteria of what is already in the buffer
    if (this->keep_collecting())
    {
//      log_os << "Started @ " << this->block_start << std::endl;
//      log_os << "keep collecting - @ " << si  << std::endl;
//      log_os << "var count " << var_count << std::endl;
//      this->write_out_buffer();
        return false;
    }

    // We are not collecting further, check if we want to assemble what is in the buffer
    if (this->do_assemble())
    {
        // if we decide to assemble, generate contig-space; modify gVCF records and buffer accordingly
        make_record();
//      log_os << "Assembling " << this->block_start << " - " << this->block_end << std::endl;
    }
    this->notify_consumer();
    this->clear();
    return true;
}

// makes the phased VCF record from the buffered sites list
void
assembly_streamer::make_record()
{
    this->construct_reference();
    this->collect_read_evidence();
    this->assemble();
    this->create_contig_records();
}

void
assembly_streamer::construct_reference()
{
    this->reference = "";
    // TODO: the following will need to be revised to handle indels!
    for (unsigned i=0; i<(_site_buffer.size()); i++)
        this->reference += _site_buffer.at(i).ref;
}

void
assembly_streamer::assemble()
{

  // THIS IS A STUB ASSEMBLER: IT ASSUMES THAT THE PREDICTOR CAN PROVIDE THE DESIRED ASSEMBLY

  int regionCount(0);
  bool inRegion(false);
  unsigned aPosInRegion=0;
  for(unsigned i=block_start;i<=block_end;++i){
    if (this->myPredictor.rt.isPayloadInRegion(i)){
      if(!inRegion){
        ++regionCount;
        if(regionCount==1){
          aPosInRegion=i;
        }
        inRegion=true;
      }
    } else {
      inRegion=false;
    }
  }
  assert(regionCount==1);

  boost::optional<std::vector<std::string> > ctgs(myPredictor.rt.isPayloadInRegion(aPosInRegion));
  if (ctgs){
    asm_contigs = *ctgs;
  } else {
    asm_contigs.clear();
  }
}

void
assembly_streamer::rescore(std::stringstream &AD)
{
    // TODO: some sort of read realignment-based scoring
    // TODO SOON: add score information here


    // add information from the alleles selected for output
    AD << 10;
    bool hasAlt(false);
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
    site_info& base = (this->_site_buffer.at(0));
    this->clear_buffer();



#if 0
    // Dummy determine two allele with most counts in observation counter
    typedef std::pair<std::string, int> allele_count_t;
    std::array<allele_count_t, 2> max_alleles = {allele_count_t("N",0),allele_count_t("N",0)};
    for (auto& obs : observations)
    {
        if (obs.second>max_alleles[0].second)
        {
                max_alleles[1]  = max_alleles[0];
                max_alleles[0]  = obs;
        }
        if (obs.second>max_alleles[1].second && max_alleles[0].first!=obs.first)
        {
                max_alleles[1] = obs;
        }
    }


    // add information from the alleles selected for output
    std::stringstream AD,alt;
    AD << this->observations[this->reference];
    bool hasAlt(false);
    for (const auto& val : max_alleles)
    {
        if (val.first==reference) continue;
        hasAlt = true;
        if (! alt.str().empty()) alt << ',';
        alt << val.first;
        AD << ',' << val.second;
    }
    if ( ! hasAlt )
    {
        alt << ".";
    }


#else
    // add information from the alleles selected for output
    std::stringstream alt;
    bool hasAlt(false);
    for (const auto& val : asm_contigs)
    {
        if (val==reference) continue;
        hasAlt = true;
        if (! alt.str().empty()) alt << ',';
        alt << val;
    }
    if ( ! hasAlt )
    {
        alt << ".";
    }

    std::stringstream AD;
    rescore(AD);

#endif
//    log_os << "max_1 " << max_alleles[0].first << "=" << max_alleles[0].second << "\n";
//    log_os << "max_2 " << max_alleles[1].first << "=" << max_alleles[1].second << "\n";

    base.phased_ref = this->reference;
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
    base.phased_alt = alt.str();
    base.phased_AD  = AD.str();

    //set any filters
    base.smod.filters.reset();
    base.smod.assemblyReason = this->myPredictor.assemblyReason;

    // specify that we are covering multiple records, and make the needed modification in the output
    base.smod.is_phased_region = true;
    base.smod.is_assembled_contig = true;

    // set more vcf record fields
    int reads_ignored   = 10;
    base.n_used_calls   = 15;
    base.n_unused_calls = 20; // second term mark all alleles that we didnt use as unused reads


    // Add in assembled record(s)
    _site_buffer.push_back(base);
}

void
assembly_streamer::
collect_read_segment_evidence(
    const read_segment& rseg)
{
    //check if we are covering the block range
//    if (static_cast<unsigned>(abs(rseg.buffer_pos-this->block_end))>rseg.read_size())
//        return;
//
//    const bam_seq bseq(rseg.get_bam_read());
//
//    // read quality checks
//    if (static_cast<int>(rseg.map_qual())<this->opt.min_single_align_score || rseg.is_invalid_realignment || !rseg.is_valid())
//    {
////                this->total_reads_unused++; // do not count filtered reads in DPF
//        return;
//    }
//
//    const int sub_start((this->block_start-rseg.buffer_pos));
//    const int sub_end((this->block_end-rseg.buffer_pos));
////#ifdef DEBUG_ASSEMBLE
////    int pad(0); // add context to extracted alleles for debugging
////    sub_start -= pad;
////    sub_end += pad;
////#endif
//    using namespace BAM_BASE;
//
//    if (sub_start>=0 && sub_end<max_read_len)
//    {
//        bool do_include(true);
//        //instead make bit array counting the number of times het pos co-occur
//        std::string sub_str("");
//        for (int t=sub_start; t<(sub_end+1); t++) //pull out substring of read
//        {
//            if (bseq.get_char(t)=='N'|| static_cast<int>(rseg.qual()[t]<this->opt.min_qscore))
//            {
//                do_include = false; // do qual check of individual bases here, kick out entire read if we dont meet cut-off
//                break;
//            }
//            sub_str+= bseq.get_char(t); //TODO use more efficient data structure here
//        }
//        if (do_include)
//        {
//            this->observations[sub_str]++;
//            total_reads++;
//        }
//        else
//            this->total_reads_unused++;
//    }
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

bool
assembly_streamer::
keep_collecting()
{
        //extend with more data structures to determine is assembly criteriea is met
        return this->myPredictor.keep_extending(this->block_start,this->block_end);
}

bool
assembly_streamer::
do_assemble()
{
        //extend with more data structures to determine is assembly criteriea is met
        return this->myPredictor.do_assemble(this->block_start,this->block_end);
}


void
assembly_streamer::clear()
{
    this->clear_buffer();
    observations.clear();
    block_start                                 = -1;
    block_end                                   = -1;
    var_count                   = 0;
    total_reads                 = 0;
    total_reads_unused          = 0;
    reference                   = "";
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
