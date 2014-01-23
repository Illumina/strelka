/*
 * codon_buffer.hh
 *
 *  Created on: Sep 13, 2013
 *      Author: mkallberg
 */

#ifndef CODON_BUFFER_HH_
#define CODON_BUFFER_HH_
#include "starling_common/starling_pos_processor_base.hh"

class Codon_buffer {
public:
    Codon_buffer();
    virtual ~Codon_buffer();
//    void set_sample(starling_pos_processor::sample_info* si);

private:
    int mysample;
//    sample_info* mysample;
};

#endif /* CODON_BUFFER_HH_ */
