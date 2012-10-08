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

/// \author Chris Saunders
///
#ifndef __STARLING_DIPLOID_INDEL_HH
#define __STARLING_DIPLOID_INDEL_HH


#include "blt_util/qscore.hh"

#include "boost/utility.hpp"


namespace STAR_DIINDEL {
    enum index_t {
        NOINDEL,
        HOM,
        HET,
        SIZE
    };

    inline
    const char*
    label(const unsigned idx){
        switch(idx){
        case NOINDEL: return "ref";
        case HOM: return "hom";
        case HET: return "het";
        default: return "xxx";
        }
    }

#if 0
    // states are the number of copies of I,R,NR
    // I=called indel allele
    // R=reference (no indels which interfere with I)
    // NR=first indel allele interfering with I
    enum index_t {
        I2,
        I1R1,
        I1NR1,
        R2,
        R1NR1,
        NR2,
        SIZE
    };

    inline
    const char*
    label(const unsigned idx){
        switch(idx){
        case I2:       return "hom";
        case I1R1:     return "het";
        case I1NR1:    return "nonref_het";
        case R2:       return "ref";
        case R1NR1:
        case NR2:      return "other";
        default:       return "xxx";
        }
    }
#endif
}



struct starling_diploid_indel : private boost::noncopyable {

    starling_diploid_indel()
        : is_indel(false), max_gt(0) {
        static const double p(1./static_cast<double>(STAR_DIINDEL::SIZE));
        static const int qp(error_prob_to_qphred((1.-p)));
        indel_qphred=qp;
        max_gt_qphred=qp;
        for(unsigned i(0);i<STAR_DIINDEL::SIZE;++i) pprob[i] = p;
    }

    bool is_indel;
    unsigned max_gt;
    int indel_qphred;
    int max_gt_qphred;
    double pprob[STAR_DIINDEL::SIZE];
};


#endif
