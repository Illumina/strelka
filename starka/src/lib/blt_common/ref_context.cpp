// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
//

/// \file
///
/// \author Chris Saunders
///


#include "blt_common/ref_context.hh"
#include "blt_util/seq_util.hh"



unsigned
get_snp_hpol_size(const pos_t pos,
                  const reference_contig_segment& ref){
    
    // count upstream repeats:
    bool is_up_repeat(false);
    char up_repeat('N');
    unsigned up_size(0);
    for(pos_t i(pos-1);i>=0;i--){
        if(is_up_repeat) {
            if(up_repeat != ref.get_base(i)) break;
        } else {
            up_repeat=ref.get_base(i);
            is_up_repeat=true;
            if(up_repeat=='N') break;
        }
        up_size++;
    }

    // count downstream repeats:
    bool is_dn_repeat(false);
    char dn_repeat('N');
    unsigned dn_size(0);
    const pos_t rs(ref.end());
    for(pos_t i(pos+1);i<rs;i++){
        if(is_dn_repeat) {
            if(dn_repeat != ref.get_base(i)) break;
        } else {
            dn_repeat=ref.get_base(i);
            is_dn_repeat=true;
            if(dn_repeat=='N') break;
        }
        dn_size++;
    }

    return 1+((dn_repeat==up_repeat) ? up_size+dn_size : std::max(up_size,dn_size) );
}



// helper class for finding interupted homopolymers:
//
struct ihpol_data {

    ihpol_data()
      : r1('N'), r2('N')
      , nr1(0), nr2(0)
      , size(0)
    {}

    // returns false when limit is reached:
    bool
    add_base(const char b) {
        if       (nr1==0) {
            r1=b;
            nr1=1;
        } else if(r1==b) {
            if((nr2>1) || (r1=='N')) return false;
            nr1++;
        } else if(nr2==0) {
            r2=b;
            nr2=1;
        } else if(r2==b) {
            if((nr1>1) || (r2=='N')) return false;
            nr2++;
        } else {
            return false;
        }
        size++;
        return true;
    }

    unsigned
    max_allele_size() const {
        return std::max(nr1,nr2);
    }

    char r1,r2;
    unsigned nr1,nr2;
    unsigned size;
};



unsigned
get_interupted_hpol_size(const pos_t pos,
                         const reference_contig_segment& ref){

    // count current base + upstream repeats:
    ihpol_data up_ihd;
    for(pos_t i(pos);i>=0;i--){
        const char base(ref.get_base(i));
        if(! up_ihd.add_base(base)) break;
    }
    // then see how far we can extend downstream:
    const pos_t rs(ref.end());
    for(pos_t i(pos+1);i<rs;i++){
        const char base(ref.get_base(i));
        if(! up_ihd.add_base(base)) break;
    }

    // count current base + downstream repeats:
    ihpol_data dn_ihd;
    for(pos_t i(pos);i<rs;i++){
        const char base(ref.get_base(i));
        if(! dn_ihd.add_base(base)) break;
    }
    // then see how far we can extend upstream:
    for(pos_t i(pos-1);i>=0;i--){
        const char base(ref.get_base(i));
        if(! dn_ihd.add_base(base)) break;
    }

    // return the longer of the two:
    return std::max(up_ihd.max_allele_size(),dn_ihd.max_allele_size());
}

