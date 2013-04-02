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
#ifndef __POSITION_SOMATIC_SNV_HH
#define __POSITION_SOMATIC_SNV_HH

#include "blt_common/snp_pos_info.hh"
#include "blt_common/position_snp_call_pprob_digt.hh"
#include "strelka/strelka_shared.hh"

#include "blt_util/qscore.hh"

#include <boost/utility.hpp>

#include <iosfwd>


//#define SOMATIC_DEBUG


namespace DDIGT {

enum index_t { SIZE = DIGT::SIZE*DIGT::SIZE };

inline
unsigned
get_state(const unsigned normal_gt,
          const unsigned tumor_gt) {
    return normal_gt+DIGT::SIZE*tumor_gt;
}

inline
void
get_digt_states(const unsigned dgt,
                unsigned& normal_gt,
                unsigned& tumor_gt) {
    normal_gt = (dgt%DIGT::SIZE);
    tumor_gt = (dgt/DIGT::SIZE);
}
}


std::ostream& operator<<(std::ostream& os,const DDIGT::index_t dgt);



struct somatic_snv_genotype : private boost::noncopyable {

    somatic_snv_genotype()
        : is_snv(false), tier(0), ref_gt(0) {}

    struct result_set {

        result_set()
            : max_gt(0), pprob(DDIGT::SIZE)
        {
            static const blt_float_t p(1./static_cast<blt_float_t>(DDIGT::SIZE));
            static const int qp(error_prob_to_qphred((1.-p)));
            snv_qphred=qp;
            snv_from_ref_qphred=qp;
            snv_from_het_qphred=qp;
            snv_from_het_loh_qphred=qp;
            snv_from_het_nonloh_qphred=qp;
            snv_from_hom_qphred=qp;
            snv_from_anyhom_qphred=qp;
            max_gt_qphred=qp;
            for(unsigned i(0); i<DDIGT::SIZE; ++i) {
                pprob[i] = p;
            }
        }

        // TODO: add marginal normal/tumor genotypes

        unsigned max_gt;
        int snv_qphred;
        int snv_from_ref_qphred;
        int snv_from_het_qphred;
        int snv_from_het_loh_qphred;
        int snv_from_het_nonloh_qphred;
        int snv_from_hom_qphred;
        int snv_from_anyhom_qphred;
        int max_gt_qphred;
        std::vector<blt_float_t> pprob;
    };

    bool is_snv;
    unsigned tier;
    unsigned ref_gt;
    result_set genome;
    result_set poly;
};



// snv call output:
//
void
write_somatic_snv_genotype(const strelka_options& opt,
                           const somatic_snv_genotype& sgt,
                           const snp_pos_info& normal_pi,
                           const snp_pos_info& tumor_pi,
                           std::ostream& os);


// object used to pre-compute priors:
struct somatic_snv_caller {

    somatic_snv_caller(const strelka_options& opt,
                       const pprob_digt_caller& pd_caller);

    //
    void
    position_somatic_snv_call(const extended_pos_info& normal_epi,
                              const extended_pos_info& tumor_epi,
                              const extended_pos_info* normal_epi_t2_ptr,
                              const extended_pos_info* tumor_epi_t2_ptr,
                              somatic_snv_genotype& sgt) const;

private:
    const blt_float_t*
    lnprior_genomic(const unsigned ref_id) const {
        return _lnprior[ref_id].genome;
    }

    const blt_float_t*
    lnprior_polymorphic(const unsigned ref_id) const {
        return _lnprior[ref_id].poly;
    }

    struct prior_set {
        blt_float_t genome[DDIGT::SIZE];
        blt_float_t poly[DDIGT::SIZE];
    };

    const strelka_options& _opt;
    prior_set _lnprior[N_BASE+1];
};

#endif
