// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#pragma once


#include "starling_common/starling_pos_processor_base.hh"

#include "boost/foreach.hpp"


#if 0
struct null_snp_pos_info {
    null_snp_pos_info() {
        for (unsigned i(0); i<(N_BASE+1); ++i) {
            pi[i].ref_base = id_to_base(i);
        }
    }

    snp_pos_info pi[N_BASE+1];
};
#endif


struct sample_pos_data {

    sample_pos_data(snp_pos_info* pi_ptr,
                    extra_position_data& epd_init,
                    const char ref_base,
                    const blt_options& opt,
                    dependent_prob_cache& dpcache,
                    const bool is_dependent_eprob,
                    const bool is_include_tier2)
        : pi((NULL==pi_ptr) ? _null_pi : *pi_ptr),
          epd(epd_init),
          n_calls(0) {

        pi.ref_base=ref_base;

        // for all but coverage-tests, we use a high-quality subset of the basecalls:
        //
        epd.good_pi.clear();
        epd.good_pi.ref_base = pi.ref_base;

        n_calls += pi.calls.size();
        BOOST_FOREACH(const base_call& bc, pi.calls) {
            if (bc.is_call_filter) {
                if (! (is_include_tier2 &&
                       bc.is_tier_specific_call_filter)) {
                    continue;
                }
            }
            epd.good_pi.calls.push_back(bc);
        }

        if (is_include_tier2) {
            n_calls += pi.tier2_calls.size();
            BOOST_FOREACH(const base_call& bc, pi.tier2_calls) {
                if (bc.is_call_filter) continue;
                epd.good_pi.calls.push_back(bc);
            }
        }

        adjust_joint_eprob(opt,dpcache,epd.good_pi,epd.dependent_eprob,is_dependent_eprob);

        n_used_calls=(epd.good_pi.calls.size());
        n_unused_calls=(n_calls-n_used_calls);
    }

private:
    //static null_snp_pos_info;
    snp_pos_info _null_pi;
public:
    snp_pos_info& pi;
    extra_position_data& epd;
    unsigned n_calls;
    unsigned n_used_calls;
    unsigned n_unused_calls;
};


struct extended_pos_data : public sample_pos_data {

    typedef sample_pos_data base_t;

    extended_pos_data(snp_pos_info* pi_ptr,
                      extra_position_data& epd_init,
                      const char ref_base,
                      const blt_options& opt,
                      dependent_prob_cache& dpc,
                      const bool is_dependent_eprob,
                      const bool is_include_tier2)
        : base_t(pi_ptr,epd_init,ref_base,opt,dpc,is_dependent_eprob,is_include_tier2)
        , good_epi(epd.good_pi,epd.dependent_eprob) {}

    extended_pos_info good_epi;
};
