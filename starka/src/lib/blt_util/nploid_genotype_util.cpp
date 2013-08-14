// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

/// \file

/// \author Chris Saunders
///
#include "blt_util/nploid_genotype_util.hh"

#include <boost/math/special_functions/binomial.hpp>



unsigned
nploid_gtype_size(const unsigned ploidy) {

    static const unsigned NB1(N_BASE-1);
    const double r(boost::math::binomial_coefficient<double>(NB1+ploidy,NB1));
    assert(r>=0);
    return static_cast<unsigned>(r);
}



nploid_info::
nploid_info(const unsigned init_ploidy)
    : _ploidy(init_ploidy), _echunk(1./static_cast<double>(init_ploidy)) {

    gtype_info gi;
    gi.label.resize(_ploidy);
    ginfo_init(_ploidy,0,0,gi);

    unsigned init_index(0);
    unsigned* init_ref(_ref_gtype);
    ref_gtype_init(_ploidy,0,0,true,init_ref,init_index);
}



// recursively derive labels and base frequencies for each genotype
void
nploid_info::
ginfo_init(const unsigned init_ploidy,
           const unsigned pli,
           const unsigned init_i,
           gtype_info& gi) {

    for (unsigned i(init_i); i<N_BASE; ++i) {
        gi.efreq_levels[i] += 1;
        gi.label[pli]=id_to_base(i);
        if ((pli+1)<init_ploidy) {
            ginfo_init(init_ploidy,pli+1,i,gi);
        } else {
            _ginfo.push_back(gi);
        }
        gi.efreq_levels[i] -= 1;
    }
}



// build an index of which genotype corresponds to the reference state,
// for each possible reference base:
void
nploid_info::
ref_gtype_init(const unsigned init_ploidy,
               const unsigned pli,
               const unsigned init_i,
               const bool is_hom,
               unsigned*& ref,
               unsigned& index) {

    for (unsigned i(init_i); i<N_BASE; ++i) {
        const bool is_iter_hom(pli==0 || (is_hom && (i==init_i)));
        if ((pli+1)<init_ploidy) {
            ref_gtype_init(init_ploidy,pli+1,i,is_iter_hom,ref,index);
        } else {
            if (is_iter_hom) *ref++ = index;
            ++index;
        }
    }
}

