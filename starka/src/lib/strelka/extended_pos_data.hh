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


#ifndef EXTENDED_POS_DATA_
#define EXTENDED_POS_DATA_


#include "starling_common/starling_pos_processor_base.hh"



#if 0
struct null_snp_pos_info {
    null_snp_pos_info() {
        for(unsigned i(0);i<(N_BASE+1);++i){
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

        const unsigned n_tier1_calls(pi.calls.size());
        for(unsigned i(0);i<n_tier1_calls;++i){
            n_calls++;
            if(pi.calls[i].is_call_filter) {
                if(not (is_include_tier2 and
                        pi.calls[i].is_tier_specific_call_filter)) {
                    continue;
                }
            }
            epd.good_pi.calls.push_back(pi.calls[i]);
        }

        if(is_include_tier2) {
            const unsigned n_tier2_calls(pi.tier2_calls.size());
            for(unsigned i(0);i<n_tier2_calls;++i){
                n_calls++;
                if(pi.tier2_calls[i].is_call_filter) continue;
                epd.good_pi.calls.push_back(pi.tier2_calls[i]);
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

#endif
