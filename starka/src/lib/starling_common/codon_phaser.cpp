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
    range       = 3;                  // phasing range we are considering
    is_in_block = false;
    het_count   = 0;
    read_len    = 100;
}

Codon_phaser::~Codon_phaser() {
    // TODO Auto-generated destructor stub
}



// Add a indel site to the phasing buffer
bool
Codon_phaser::add_site(int i) {
    //currently for indel sites we
}
// Add a SNP site to the phasing buffer
bool
Codon_phaser::add_site(site_info& si) {
    // case: extending block with het call, update block_end position
    if (si.is_het()) {
        if (!is_in_block)
            block_start = si.pos;
        is_in_block = true;
        block_end = si.pos;
        buffer.push_back(si);
        het_count ++;
        return false;
    }

    // case: extending block with none-het call based on the phasing range
    if (is_in_block && (si.pos-block_end+1)<this->range) {
        buffer.push_back(si);
        return false;
    }

    // case: setup the phased record and write out
    if (het_count>1) {
        log_os << "!!!het count " << het_count << "\n";
        this->write_out_buffer();
        make_record();
    }
    return true;
}

// makes the phased VCF record from the
// buffered sites list
site_info
Codon_phaser::make_record() {
    site_info base = buffer.at(0);
    log_os << "!!!Im in make_record - block_end " << (block_end+1) << "\n";
    int buffer_start = (block_start-this->read_len);
    int buffer_end = (block_start);
    int total_reads(0);

    // extract evidence for all reads that span the entire phasing range
    for (int i=buffer_start;i<buffer_end;i++){
        read_segment_iter ri(read_buffer->get_pos_read_segment_iter(i));
        read_segment_iter::ret_val r;
        while (true) {
            r=ri.get_ptr();
            if (NULL==r.first) break;
            const read_segment& rseg(r.first->get_segment(r.second));
//            log_os << "read qual " << static_cast<int>(rseg.qual()[10]) << "\n";
            const bam_seq bseq(rseg.get_bam_read());

            //do mapq check for read
            //rseg.map_qual();

            int sub_start((this->block_start-rseg.buffer_pos));
            int sub_end((this->block_end-rseg.buffer_pos));

            #ifdef DEBUG_CODON
                int pad(1); // add
                sub_start -= pad;
                sub_end += pad;
            #endif

            if (sub_start>0){
                std::string sub_str("");

                for (int t=sub_start;t<sub_end+1;t++)
                    sub_str+= bseq.get_char(t);
                    // do qual check of individual bases

                #ifdef DEBUG_CODON
                    log_os << "substart " << sub_start << "\n";
                    log_os << "subend " << sub_end << "\n";
                    log_os << "substr " << sub_str << "\n";
                    log_os << "read key " << rseg.key() << "\n";
                    log_os << "read pos " << rseg.buffer_pos << "\n";
                    log_os << "read seq " << bseq << "\n\n";
                #endif
                total_reads++;
            }
            ri.next();
        }
    }
    log_os << "total reads " << total_reads << "\n";

//    for (unsigned i=0; i<buffer.size(); i++) {
//        call = buffer.at(i);
//        log_os << call << "\n";
//        log_os << call.ref << "\n";
//
//        if (i==0)
//            base = call;
//        else {
//        }
//        if (call.pos==block_end)
//            break;
//    }
    return base;
}

void Codon_phaser::collect_read_evidence(){

}

void
Codon_phaser::clear_buffer() {
    buffer.clear();
    this->observations.clear();
    block_end       = -1;
    is_in_block     = false;
    het_count       = 0;
}

void
Codon_phaser::clear_read_buffer(int i){
    // clear read buffer up to next feasible position
}


void
Codon_phaser::write_out_buffer() {
    for (std::vector<site_info>::iterator it = buffer.begin(); it != buffer.end(); ++it) {
        log_os << *it << "\n";
    }
}

void
Codon_phaser::write_out_buffer() {
    for (std::vector<site_info>::iterator it = buffer.begin(); it != buffer.end(); ++it) {
        log_os << *it << "\n";
    }
}
