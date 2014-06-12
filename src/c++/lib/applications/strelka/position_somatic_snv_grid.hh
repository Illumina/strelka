// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

/// \file

/// variation on the original strowman snv caller -- implements a
/// compile-time specified grid in allele frequency space and requires
/// similar frequency as definition of non-somatic.
///

/// \author Chris Saunders
///

#pragma once

#if 0

#include "blt_common/snp_pos_info.hh"
#include "blt_common/position_snp_call_pprob_digt.hh"
#include "strelka/strelka_shared.hh"

#include "position_somatic_snv_grid_shared.hh"

#include "blt_util/qscore.hh"

#include <boost/utility.hpp>

#include <iosfwd>


//#define SOMATIC_DEBUG


namespace DIGT_GRID {

// HET_RES is the number of points sampled between 0 and 0.5 on
// the continuous frequency scale. Thus a fully sampled axis will
// be sampled HET_RES*2+3 times.
//
enum constants { HET_RES = 4,
                 HET_COUNT = HET_RES*2+1,
                 HET_SIZE = DIGT::SIZE-N_BASE
               };

enum index_t { SIZE = N_BASE+HET_SIZE*HET_COUNT };


inline
unsigned
get_digt_state (const unsigned state) {
    if (state<N_BASE) return state;
    return N_BASE+((state-N_BASE)%HET_SIZE);
}

inline
unsigned
get_het_count(const unsigned state) {
    if (state<N_BASE) return 0;
    return (state-N_BASE)/HET_SIZE;
}
}


namespace DDIGT_GRID {

enum index_t { SIZE = DIGT_GRID::SIZE*DIGT_GRID::SIZE };

inline
unsigned
get_state(const unsigned normal_gt,
          const unsigned tumor_gt) {
    return normal_gt+DIGT_GRID::SIZE*tumor_gt;
}

inline
void
get_digt_grid_states(const unsigned dgt,
                     unsigned& normal_gt,
                     unsigned& tumor_gt) {
    normal_gt = (dgt%DIGT_GRID::SIZE);
    tumor_gt = (dgt/DIGT_GRID::SIZE);
}


struct is_nonsom_maker_t {
    is_nonsom_maker_t();

    std::vector<bool> val;
};
}


std::ostream& operator<<(std::ostream& os,const DDIGT_GRID::index_t dgt);


// object used to pre-compute priors:
struct somatic_snv_caller_grid {

    somatic_snv_caller_grid(const strelka_options& opt,
                            const pprob_digt_caller& pd_caller);

    //
    void
    position_somatic_snv_call(const extended_pos_info& normal_epi,
                              const extended_pos_info& tumor_epi,
                              const extended_pos_info* normal_epi_t2_ptr,
                              const extended_pos_info* tumor_epi_t2_ptr,
                              somatic_snv_genotype_grid& sgt) const;

private:
    const std::vector<blt_float_t>&
    lnprior_genomic(const unsigned ref_id) const {
        return _lnprior[ref_id].normal;
    }

    const std::vector<blt_float_t>&
    lnprior_polymorphic(const unsigned ref_id) const {
        return _lnprior[ref_id].normal_poly;
    }

    struct prior_set {
        prior_set()
            : normal(DIGT_GRID::SIZE)
            , normal_poly(DIGT_GRID::SIZE) {}

        std::vector<blt_float_t> normal;
        std::vector<blt_float_t> normal_poly;
    };

    const strelka_options& _opt;
    const pprob_digt_caller& _pd_caller;
    prior_set _lnprior[N_BASE+1];
    blt_float_t _ln_som_match;
    blt_float_t _ln_som_mismatch;
};


// snv call output:
//
void
write_somatic_snv_genotype_grid(const strelka_options& opt,
                                const somatic_snv_genotype_grid& sgt,
                                const snp_pos_info& normal_pi,
                                const snp_pos_info& tumor_pi,
                                std::ostream& os);

#endif
