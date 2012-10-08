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
#include "blt_common/position_snp_call_bacon.hh"
#include "blt_util/seq_util.hh"

#include <cassert>
#include <cmath>



///
///
void
position_snp_call_bacon(const snp_pos_info& pi,
                        bacon_scores& bas){

    if(pi.ref_base=='N') return;

    const unsigned n_calls(pi.calls.size());

    double hom_score[N_BASE];
    for(unsigned b(0);b<N_BASE;++b) hom_score[b] = 0.;

    static const double one_third(1./3.);
    static const double log_one_third(std::log(one_third));
    for(unsigned i(0);i<n_calls;++i){
        const double eprob(pi.calls[i].error_prob());

        const double val0(std::log(1.-eprob*one_third)+log_one_third);
        const double val1(std::log(eprob)+log_one_third);

        const uint8_t obs_id(pi.calls[i].base_id);
        assert(obs_id!=BASE_ID::ANY);
        for(unsigned b(0);b<N_BASE;++b){
            if(b==obs_id){
                hom_score[b] += val1;
            } else {
                hom_score[b] += val0;
            }
        }
    }

    static const double minv_log_ten(-1./std::log(10.));
    for(unsigned b(0);b<N_BASE;++b) hom_score[b] *= minv_log_ten;

    // get max and max2:
    bas.max_base_id=0;
    bas.max2_base_id=1;
    for(unsigned b(1);b<N_BASE;++b){
        if(hom_score[b] > hom_score[bas.max_base_id]){
            bas.max2_base_id = bas.max_base_id;
            bas.max_base_id = b;
        } else if(hom_score[b] > hom_score[bas.max2_base_id]) {
            bas.max2_base_id = b;
        }
    }

    double max3(0);
    bool is_first(true);
    for(unsigned b(0);b<N_BASE;++b){
        if(b==bas.max_base_id || b==bas.max2_base_id) continue;
        if(is_first || hom_score[b] > max3){
            max3 = hom_score[b];
            is_first=false;
        }
    }
    bas.max_score = hom_score[bas.max_base_id] - max3;
    bas.max2_score = hom_score[bas.max2_base_id] - max3;

    // make special-case behavior consistent:
    if(bas.max_score == bas.max2_score){
        if(bas.max2_base_id == base_to_id(pi.ref_base)) {
            std::swap(bas.max_base_id,bas.max2_base_id);
        }
    }

    bas.is_valid = true;
}
