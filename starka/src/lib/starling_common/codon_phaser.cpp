/*
 * Codon_phaser.cpp
 *
 *  Created on: Aug 10, 2013
 *  Author: Morten Kallberg
 */

#include "codon_phaser.hh"
#include <vector>

#define DEBUG_CODON


#ifdef DEBUG_CODON
#include "blt_util/log.hh"
#endif

Codon_phaser::Codon_phaser() {
    block_end = -1;
    range = 3;                  // phasing range we are considering
    is_in_block = false;
    het_count = 0;
}

Codon_phaser::~Codon_phaser() {
    // TODO Auto-generated destructor stub
}



// Add a site to the phasing buffer
bool
Codon_phaser::add_site(site_info& si) {
    // case: extending block with het call, update block_end position
    if (si.is_het()) {
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
    site_info base; //call to me modified
    site_info call;
    log_os << "!!!Im in make_record - block_end " << (block_end+1) << "\n";
//    read_segment_iter my_iter = read_buffer->get_pos_read_segment_iter(234276351);
    for (unsigned i=0; i<buffer.size(); i++) {
        call = buffer.at(i);
        log_os << call << "\n";
        log_os << call.ref << "\n";

        if (i==0)
            base = call;
        else {
        }
        if (call.pos==block_end)
            break;
    }
    return base;
}

//void
//starling_pos_processor_base::
//buffer_codon_reads(const pos_t pos) {
//    static const bool is_include_submapped(false);
//
//    for (unsigned s(0); s<_n_samples; ++s) {
//        read_segment_iter ri(sample(s).read_buff.get_pos_read_segment_iter(pos));
//        read_segment_iter::ret_val r;
//        while (true) {
//            r=ri.get_ptr();
//            if (NULL==r.first) break;
//            const read_segment& rseg(r.first->get_segment(r.second));
//            if (is_include_submapped || rseg.is_treated_as_anytier_mapping()) {
//                pileup_read_segment(rseg,s);
//            }
//            ri.next();
//        }
//    }
//}


void
Codon_phaser::write_out_buffer() {

    for (std::vector<site_info>::iterator it = buffer.begin(); it != buffer.end(); ++it) {
        log_os << *it << "\n";
    }
}

void
Codon_phaser::clear_buffer() {
    buffer.clear();
    block_end = -1;
    is_in_block = false;
    het_count = 0;
}
