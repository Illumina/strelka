/*
 * Codon_phaser.cpp
 *
 *  Created on: Sep 10, 2013
 *  Author: Morten Kallberg
 */

#include "codon_phaser.hh"
#include <vector>

#define DEBUG_CODON


#ifdef DEBUG_CODON
#include "blt_util/log.hh"
#endif

Codon_phaser::Codon_phaser() {
    block_start   = -1;
    block_end   = -1;
    range       = 3;                  // phasing range we are considering, may be determined dynamically
    is_in_block = false;
    het_count   = 0;
    read_len    = 100;                // this options need set externally
    min_mapq    = 20;
    min_baseq   = 17;
    phase_indels= false;              // if false we break the block when encountering an indel
    this->clear_buffer();
}

//Codon_phaser::~Codon_phaser() {}

// Add a indel site to the phasing buffer
//bool
//Codon_phaser::add_site_indel(int i) {
//    //currently for indel sites we
//}

// Add a SNP site to the phasing buffer
bool
Codon_phaser::add_site(site_info& si) {
    // case: extending block with het call, update block_end position
    buffer.push_back(si);

    if (si.is_het()) {
        if (!is_in_block)
            block_start = si.pos;
        is_in_block = true;
        block_end = si.pos;
        het_count ++;
        #ifdef DEBUG_CODON
//            log_os << "starting block @ " << (this->block_start+1) << " with " << si << "\n";
        #endif
        return false;
    }

    // case: extending block with none-het call based on the phasing range
    if (is_in_block && (si.pos-block_end+1)<this->range) {
        #ifdef DEBUG_CODON
//            log_os << "Extending block with @ " << (this->block_start+1) << " with " << si << "\n";
        #endif
        return false;
    }

    // case: setup the phased record and write out
    if (het_count>1) {
//        log_os << "!!!het count " << het_count << "\n";
        this->write_out_buffer();
        make_record();
    }
    return true;
}

void
Codon_phaser::construct_reference(){
    this->reference = "";
    for (unsigned i=0;i<this->buffer.size()-(this->range-1);i++)
        this->reference += buffer.at(i).ref;
}

void
Codon_phaser::create_phased_record() {
    if (this->total_reads<10){  // some initial minimum conditions, look for at least 10 spanning reads support
        // set flag saying to little evidence to phase
        return;
    }

    //Decide if we accept the novel alleles, very hacky criteria for now
    //at least 90% of the reads have to support a diploid model
    //TODO unphased genotype corresponds to as many phased genotypes as there are permutations of
    //the observed alleles. Thus, for a given unphased genotyping G1, . . . ,Gn,
    //we need to to calculate probability of sampling a particular unphased genotype combinations
    //given a set of allele frequencies: P(G1, . . . ,Gn|f1, . . . , fk)



    site_info &base = (this->buffer.at(0));
    base.region_ref = this->reference;
    base.smod.is_phased_region = true;
    base.n_used_calls = this->total_reads;
    base.n_unused_calls = this->total_reads_unused;

    // case we want to report the phased record instead
    this->buffer.clear();
    this->buffer.push_back(base);
}

// makes the phased VCF record from the buffered sites list
void
Codon_phaser::make_record() {
    this->construct_reference();
    this->collect_read_evidence();
    this->create_phased_record();

    #ifdef DEBUG_CODON
        this->write_out_alleles();
    #endif
}

void Codon_phaser::collect_read_evidence(){
    int buffer_start = (block_start-this->read_len);
    int buffer_end = (block_start);
    // extract evidence for all reads that span the entire phasing range
    for (int i=buffer_start;i<buffer_end;i++){
        read_segment_iter ri(read_buffer->get_pos_read_segment_iter(i));
        read_segment_iter::ret_val r;
        while (true) {
            bool do_include(true);
            r=ri.get_ptr();
            if (NULL==r.first) break;
            const read_segment& rseg(r.first->get_segment(r.second));
            const bam_seq bseq(rseg.get_bam_read());

            // read quality checks
            if (rseg.is_invalid_realignment || !rseg.is_valid()) break;
            if (static_cast<int>(rseg.map_qual())<this->min_mapq) break;

            int sub_start((this->block_start-rseg.buffer_pos));
            int sub_end((this->block_end-rseg.buffer_pos));
            #ifdef DEBUG_CODON
                int pad(0); // add context to extracted alleles for debugging
                sub_start -= pad;
                sub_end += pad;
            #endif
            using namespace BAM_BASE;

            if (sub_start>=0 && sub_end<read_len){
                //instead make bit array counting the number of times het pos co-occur
                //explore the space of all phas-possibilities?
                std::string sub_str("");
                for (int t=sub_start;t<(sub_end+1);t++){ //pull out substring of read
                    if (bseq.get_char(t)=='N'|| static_cast<int>(rseg.qual()[t]<this->min_baseq)){
                        do_include = false; // do qual check of individual bases here, kick out entire read if we dont meet cut-off
                        break;
                    }
                    sub_str+= bseq.get_char(t); //TODO use more efficient data structure here
                }
                #ifdef DEBUG_CODON
//                    log_os << "substart " << sub_start << "\n";
//                    log_os << "subend " << sub_end << "\n";
//                    log_os << "substr " << sub_str << "\n";
//                    log_os << "read key " << rseg.key() << "\n";
//                    log_os << "read pos " << rseg.buffer_pos << "\n";
//                    log_os << "do_include " << do_include << "\n";
//                    log_os << "read seq " << bseq << "\n\n";
                #endif
                if(do_include){
                    this->observations[sub_str]++;

                    // maintain max keys
                    if (this->max_alleles[0]==""){
                        this->max_alleles[0] = sub_str;
                        this->max_alleles[1] = sub_str;
                    }
                    else if (this->observations[sub_str]> this->observations[this->max_alleles[0]]){
                        this->max_alleles[1] = this->max_alleles[0];
                        this->max_alleles[0] = sub_str;
                    }
                    else if(this->observations[sub_str]> this->observations[this->max_alleles[1]])
                        this->max_alleles[1] = sub_str;

                    total_reads++;
                }
                else
                    this->total_reads_unused++;
            }
            ri.next();
        }
    }

    #ifdef DEBUG_CODON
        log_os << "total reads " << total_reads << "\n";
    #endif
}

void
Codon_phaser::clear_buffer() {
    buffer.clear();
    this->observations.clear();
    block_end               = -1;
    is_in_block             = false;
    het_count               = 0;
    this->total_reads       = 0;
    this->total_reads_unused = 0;
    this->reference          = "";
    for (unsigned i(0);i<2;i++)
        this->max_alleles[i] = "";
}

void
Codon_phaser::clear_read_buffer(int i){
    // clear read buffer up to next feasible position
    log_os << i << "\n";
}


void
Codon_phaser::write_out_buffer() {
    for (std::vector<site_info>::iterator it = buffer.begin(); it != buffer.end(); ++it)
        log_os << *it << "\n";
}

void
Codon_phaser::write_out_alleles() {
    log_os << "Ref: " << this->reference << "=" << this->observations[this->reference] << "\n";
    log_os << "Alt: ";
    for (allele_map::const_iterator it = this->observations.begin(); it != this->observations.end(); ++it) { // only normalize the features that are needed
        if (it->first!=this->reference)
            log_os << it->first << "=" << this->observations[it->first] << " ";
    }
    log_os << "\n";
}
