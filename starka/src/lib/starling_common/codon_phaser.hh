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
    bool add_site(site_info& si);       // add site to buffer
    void clear_buffer(void);            // clear site buffer
    site_info make_record(void);        // make phased record
    void write_out_buffer();            // debugging feature, print current buffer to std
    void clear_read_buffer(int pos);    // free up read that are no longer in phasing evidence, up to and including this position
    void collect_read_evidence();       // fill in allele counter
    bool is_in_block;                   // Are we currently in a phasing block
    std::vector<site_info> buffer;      // buffer of het snp calls
    starling_read_buffer *read_buffer;  // pass along the relevant read-buffer
private:
    int block_start;                    // position of first added het site to block
    int block_end;                      // position of last added het site to block
    int range;                          // phasing window considered
    int het_count;                      // total hets observed in buffer
    int read_len;                       // the length of the input reads
    int previous_clear;                 // cleared buffer up to this site
    std::string reference;              // the phased allele reference
    std::map<std::string,int> observations;
};
#endif /* CODONPHASER_HH_ */
