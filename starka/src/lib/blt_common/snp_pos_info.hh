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
#ifndef __SNP_POS_INFO_HH
#define __SNP_POS_INFO_HH

#include "blt_common/hapscore.hh"
#include "blt_common/ranksum.hh"
#include "blt_util/qscore.hh"
#include "blt_util/seq_util.hh"

#include "stdint.h"

#include <iosfwd>
#include <vector>

//#define BC_DEBUG



struct base_call {
    base_call(const uint8_t init_base_id,
              const uint8_t init_qscore,
              const bool init_ifs,
#ifdef BC_DEBUG
              const uint16_t init_read_pos,
              const uint16_t init_read_size,
#else
              const uint16_t,
              const uint16_t,
#endif
              const bool init_is_call_filter,
              const bool init_is_neighbor_mismatch,
              const bool init_is_tscf=false
             )
        : qscore(init_qscore),
#ifdef BC_DEBUG
          read_pos(init_read_pos), read_size(init_read_size),
#endif
          base_id(init_base_id),
          is_fwd_strand(init_ifs),
          is_neighbor_mismatch(init_is_neighbor_mismatch),
          is_call_filter(init_is_call_filter),
          is_tier_specific_call_filter(init_is_tscf)
    {
        qphred_cache::qscore_check(qscore,"basecall quality");
    }

    // pull quality value transformations from caching functions:
    double
    error_prob() const {
        return qphred_to_error_prob(static_cast<int>(qscore));
    }

    double
    ln_error_prob() const {
        return qphred_to_ln_error_prob(static_cast<int>(qscore));
    }

    double
    ln_comp_error_prob() const {
        return qphred_to_ln_comp_error_prob(static_cast<int>(qscore));
    }

    uint16_t
    get_qscore() const { return qscore; }

private:
    uint16_t qscore:8;
public:
    uint16_t base_id:4;
    uint16_t is_fwd_strand:1;
    uint16_t is_neighbor_mismatch:1;
    uint16_t is_call_filter:1; // filtered from snp-calling
    uint16_t is_tier_specific_call_filter:1;
#ifdef BC_DEBUG
    uint16_t read_pos; // zero-indexed cycle number
    uint16_t read_size;
#endif
};

std::ostream& operator<<(std::ostream& os,const base_call& bc);



struct snp_pos_info {

    snp_pos_info() { clear(); }

    void
    clear() {
        ref_base='N';
        calls.clear();
        tier2_calls.clear();
        is_n_ref_warn=false;
        n_spandel=0;
        n_mapq=0;
        cumm_mapq=0;
        mq_ranksum.clear();
        baseq_ranksum.clear();
        read_pos_ranksum.clear();
        hap_set.clear();
    }

    template <typename T>
    void
    get_known_counts(T& base_count,
                     const int min_qscore) const {

        for(unsigned i(0); i<N_BASE; ++i) base_count[i] = 0;

        const unsigned n_calls(calls.size());
        for(unsigned i(0); i<n_calls; ++i) {
            if(calls[i].base_id==BASE_ID::ANY) continue;
            if(calls[i].get_qscore()<min_qscore) continue;
            base_count[calls[i].base_id]++;
        }
    }


    void
    set_ref_base(char base){
        ref_base = base;
        mq_ranksum.set_ref_base(base);
        baseq_ranksum.set_ref_base(base);
        read_pos_ranksum.set_ref_base(base);
    }

    char
    get_ref_base(){
        return ref_base;
    }

    // returns the RMS of the read mapQs
    double
    get_rms_mq();

    // returns the read-position rank sum
    double
    get_read_pos_ranksum();

    // returns the mapQ rank sum
    double
    get_mq_ranksum();

    // return the baseq rank sum
    double
    get_baseq_ranksum();

    void
    print_known_counts(std::ostream& os,
                       const int min_qscore) const;

    void
    print_known_qscore(std::ostream& os,
                       const int min_qscore) const;

public:
    char ref_base; // always fwd-strand base
    bool is_n_ref_warn;
    std::vector<base_call> calls;
    std::vector<base_call> tier2_calls; // call not passing stringent quality criteria
    // number of spanning deletions crossing the site:
    unsigned n_spandel;
    // number of submapped reads crossing the site.
    // note this could be usable,filtered or spanning deletion,
    // all submapped reads get counted here:
    unsigned n_submapped;

    // number of mapq observations, meaning coverage at position
    unsigned n_mapq;
    // sum of mapq for all reads at this position
    int cumm_mapq;

    hap_set_t hap_set;

    //for calculating various rank-sum statistics
    ranksum mq_ranksum;
    ranksum baseq_ranksum;
    ranksum read_pos_ranksum;
};

std::ostream& operator<<(std::ostream& os,const snp_pos_info& pci);



struct extended_pos_info {
    extended_pos_info(const snp_pos_info& pi_init,
                      const std::vector<float>& de_init)
        : pi(pi_init)
        , de(de_init) {}

    const snp_pos_info& pi;
    const std::vector<float>& de;
};

#endif
