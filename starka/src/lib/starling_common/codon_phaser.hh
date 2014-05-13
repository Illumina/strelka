/*
 * Codonphaser.hh
 *
 *  Test for codon-phasing.
 *
 *  Created on: Aug 10, 2013
 *  Author: Morten Kallberg
 */

#ifndef CODONPHASER_HH_
#define CODONPHASER_HH_

#include "starling_common/gvcf_locus_info.hh"
#include "starling_common/starling_read_buffer.hh"

class Codon_phaser {
public:
    Codon_phaser();
    virtual ~Codon_phaser();

public:

    bool add_site(site_info& si);  // Add site to buffer
    void clear_buffer(void);       // clear site buffer
    site_info make_record(void);   // make phased record
    void write_out_buffer();       // debugging feature, print current buffer to std
    bool is_in_block;              // Are we currently in a phasing block
    std::vector<site_info> buffer;
    starling_read_buffer *read_buffer; // pass along the relevant read-buffers
private:
    int block_end;                 // position of last added het site to block
    int range;
    int het_count;
};

#endif /* CODONPHASER_HH_ */
