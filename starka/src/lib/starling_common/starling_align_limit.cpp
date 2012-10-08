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

#include "starling_common/starling_align_limit.hh"

#include "boost/math/special_functions/binomial.hpp"

#include <cmath>



#if 0
// Calculate the maximum number of candidate alignments for a given
// set of togglable indels and a toggle depth.
//
static
float
max_candidate_alignment_count(const unsigned n_indel,
                              const unsigned n_toggle) {

    const float n(n_indel);

    float sum(1);
    for(unsigned i(0);i<n_toggle;++i){
        const float k(i+1);
        sum += std::pow(2.,k)*boost::math::binomial_coefficient<float>(n,k);
    }
    return sum;
}
#endif



// Calculate the number of indel toggles which cannot exceed the
// maximum candidate alignment count.
//
static
unsigned
max_candidate_alignment_toggle(const unsigned n_indel,
                               const unsigned max_alignments) {

    const float max(max_alignments);

    float sum(1.);
    for(unsigned i(0);i<n_indel;++i){
        const unsigned k(i+1);
        sum += std::pow(static_cast<float>(2),static_cast<float>(k))*boost::math::binomial_coefficient<float>(n_indel,k);
        if(sum>max) return i;
   }
    return n_indel;
}



starling_align_limit::
starling_align_limit(const unsigned mac) {

    static const unsigned max_indels(100);
    for(unsigned i(0);i<max_indels;++i){
        const unsigned mt(max_candidate_alignment_toggle(i,mac));
        if((i>1) && (mt<2)) break;
        _max_toggle.push_back(max_candidate_alignment_toggle(i,mac));
    }
}
