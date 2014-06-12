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

/// variation on the original strawman snv caller -- implements a
/// compile-time specified grid in allele frequency space and requires
/// similar frequency as definition of non-somatic.
///

/// \author Chris Saunders
///

#pragma once

#include "extended_pos_data.hh"
#include "position_somatic_snv_grid_shared.hh"

#include "blt_common/snp_pos_info.hh"
#include "blt_common/position_snp_call_pprob_digt.hh"
#include "blt_util/qscore.hh"
#include "strelka/strelka_shared.hh"

#include <boost/utility.hpp>

#include <iosfwd>


//#define SOMATIC_DEBUG


namespace DIGT_SGRID {

// HET_RES is the number of points sampled between 0 and 0.5 on
// the continuous frequency scale. Thus a fully sampled axis will
// be sampled HET_RES*2+3 times.
//
// Note that the single-strand tests are only implemented on the
// half axis containing the reference allele as a major
// allele. The "STRAND_STATE_SIZE" appended to the PRESTAND_SIZE
// represents these additional strand-specific error states.
//
enum constants { HET_RES = 4,
                 HET_COUNT = HET_RES*2+1,
                 HET_SIZE = DIGT::SIZE-N_BASE,
                 STRAND_COUNT = HET_RES,
                 STRAND_SIZE = HET_SIZE/2,
                 HET_STATE_SIZE = HET_SIZE*HET_COUNT,
                 PRESTRAND_SIZE = N_BASE+HET_STATE_SIZE,
                 STRAND_STATE_SIZE = STRAND_COUNT*STRAND_SIZE
               };

enum index_t { SIZE = PRESTRAND_SIZE+STRAND_STATE_SIZE };

// generates fast lookup-table to translate between stranded ref-only
// het states and DIGT het states:
//
struct strand_state_tables {
    strand_state_tables();

    unsigned digt_state[N_BASE][STRAND_SIZE];
};

extern const strand_state_tables stables;

inline
unsigned
get_het_count(const unsigned state) {
    if (state<N_BASE)         return 0;
    if (state<PRESTRAND_SIZE) return (state-N_BASE)/HET_SIZE;
    return (state-PRESTRAND_SIZE)/STRAND_SIZE;
}

inline
bool
is_strand_state(const unsigned state) {
    return (state>=PRESTRAND_SIZE);
}

inline
unsigned
get_strand_state(const unsigned state) {
    return (state-PRESTRAND_SIZE)%STRAND_SIZE;
}

inline
unsigned
get_digt_state(const unsigned state,
               const unsigned ref_base) {
    if (state<N_BASE)         return state;
    if (state<PRESTRAND_SIZE) return N_BASE+((state-N_BASE)%HET_SIZE);
    return stables.digt_state[ref_base][get_strand_state(state)];
}


// write only most representative genotype, ie "GT", for each state
void
write_state(const DIGT_SGRID::index_t gt,
            const unsigned ref_gt,
            std::ostream& os);

// write genotype, heterozygous id, and whether this is a single-strand state
void
write_full_state(const DIGT_SGRID::index_t gt,
                 const unsigned ref_gt,
                 std::ostream& os);
}


namespace DDIGT_SGRID {

enum constants { PRESTRAND_SIZE = DIGT_SGRID::PRESTRAND_SIZE*DIGT_SGRID::PRESTRAND_SIZE };
enum index_t { SIZE = PRESTRAND_SIZE+DIGT_SGRID::STRAND_STATE_SIZE };

inline
unsigned
get_state(const unsigned normal_gt,
          const unsigned tumor_gt) {
    if (normal_gt<DIGT_SGRID::PRESTRAND_SIZE) return normal_gt+DIGT_SGRID::PRESTRAND_SIZE*tumor_gt;
    return PRESTRAND_SIZE+normal_gt-DIGT_SGRID::PRESTRAND_SIZE;
}

inline
void
get_digt_grid_states(const unsigned dgt,
                     unsigned& normal_gt,
                     unsigned& tumor_gt) {
    if (dgt<PRESTRAND_SIZE) {
        normal_gt = (dgt%DIGT_SGRID::PRESTRAND_SIZE);
        tumor_gt = (dgt/DIGT_SGRID::PRESTRAND_SIZE);
    } else {
        normal_gt= dgt+DIGT_SGRID::PRESTRAND_SIZE-PRESTRAND_SIZE;
        tumor_gt=normal_gt;
    }
}

// writes state to pattern: "AA->AG"
void
write_state(const DDIGT_SGRID::index_t dgt,
            const unsigned ref_gt,
            std::ostream& os);

// writes state to pattern: "AA_0->AG_5_strand", using a unique id for
// each heterozygous state. Writing actual allele frequencies instead of just
// ids TBD
void
write_full_state(const DDIGT_SGRID::index_t dgt,
                 const unsigned ref_gt,
                 std::ostream& os);

void
write_alt_alleles(const DDIGT_SGRID::index_t dgt,
                  const unsigned ref_gt,
                  std::ostream& os);

struct is_nonsom_maker_t {
    is_nonsom_maker_t();

    std::vector<bool> val;
};
}


std::ostream& operator<<(std::ostream& os,const DDIGT_SGRID::index_t dgt);


// object used to pre-compute priors:
struct somatic_snv_caller_strand_grid {

    somatic_snv_caller_strand_grid(const strelka_options& opt,
                                   const pprob_digt_caller& pd_caller);

    //
    void
    position_somatic_snv_call(
        const extended_pos_info& normal_epi,
        const extended_pos_info& tumor_epi,
        const extended_pos_info* normal_epi_t2_ptr,
        const extended_pos_info* tumor_epi_t2_ptr,
        const bool isComputeNonSomatic,
        somatic_snv_genotype_grid& sgt) const;

    // compute a lot of prior information for various alternate
    // versions of the method -- we don't actually need all of this for any one computation:
    //
    struct prior_set {
        prior_set()
            : normal(DIGT_SGRID::SIZE)
            , somatic_marginal(DIGT_SGRID::SIZE)
            , normal_poly(DIGT_SGRID::SIZE)
            , somatic_marginal_poly(DIGT_SGRID::SIZE)
            , normal_nostrand(DIGT_SGRID::SIZE)
            , normal_poly_nostrand(DIGT_SGRID::SIZE)
        {}

        typedef std::vector<blt_float_t> prior_t;

        prior_t normal;
        prior_t somatic_marginal;
        prior_t normal_poly;
        prior_t somatic_marginal_poly;

        // added to support somatic gVCF:
        prior_t normal_nostrand;
        prior_t normal_poly_nostrand;
    };

private:

    const prior_set&
    get_prior_set(const unsigned ref_id) const {
        return _lnprior[ref_id];
    }

    const std::vector<blt_float_t>&
    lnprior_genomic(const unsigned ref_id) const {
        return _lnprior[ref_id].normal;
    }

    const std::vector<blt_float_t>&
    lnprior_polymorphic(const unsigned ref_id) const {
        return _lnprior[ref_id].normal_poly;
    }

    const strelka_options& _opt;
    prior_set _lnprior[N_BASE+1];
    blt_float_t _ln_som_match;
    blt_float_t _ln_som_mismatch;
};



// vcf output:
//
void
write_vcf_somatic_snv_genotype_strand_grid(
    const strelka_options& opt,
    const somatic_snv_genotype_grid& sgt,
    const bool is_write_nqss,
    const extended_pos_data& n1_epd,
    const extended_pos_data& t1_epd,
    const extended_pos_data& n2_epd,
    const extended_pos_data& t2_epd,
    std::ostream& os);
