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

/// \author Chris Saunders
///

#pragma once

#include "somatic_indel_call.hh"
#include "strelka_shared.hh"

#include "starling_common/starling_indel_call_pprob_digt.hh"

#include "boost/utility.hpp"



namespace STAR_DIINDEL_GRID
{

// HET_RES is the number of frequency points sampled on a
// half-axis (that is between 0-0.5)
//
enum constants { HET_RES = 4,
                 HET_COUNT = HET_RES*2+1,
                 HET_SIZE = 1,
                 HOM_SIZE = 2
               };

// The first three states are designed to overlap with
// STAR_DIINDEL (ie. conventional diploid indel model), after this
// the grid frequency (ie. approximations of continuous frequency)
// states are added. The grid states are treated just like the
// STAR_DIINDEL het state for certain purposes (printing, for
// instance)
//
enum index_t { SIZE = HOM_SIZE+HET_SIZE*HET_COUNT };

inline
unsigned
get_star_diindel_state (const unsigned state)
{
    if (state<HOM_SIZE) return state;
    return STAR_DIINDEL::HET;
}

inline
unsigned
get_het_count(const unsigned state)
{
    if (state<HOM_SIZE) return 0;
    return (state-HOM_SIZE)/HET_SIZE;
}
}


namespace DDIINDEL_GRID
{

// This is the SIZE to represent the normal allele freq X tumor
// allele freq space which we calculate the final posteriour
// probabilty in.
//
enum index_t { SIZE = STAR_DIINDEL_GRID::SIZE*STAR_DIINDEL_GRID::SIZE };

inline
unsigned
get_state(const unsigned normal_gt,
          const unsigned tumor_gt)
{
    return normal_gt+STAR_DIINDEL_GRID::SIZE*tumor_gt;
}

inline
void
get_sdiindel_grid_states(const unsigned dgt,
                         unsigned& normal_gt,
                         unsigned& tumor_gt)
{
    normal_gt = (dgt%STAR_DIINDEL_GRID::SIZE);
    tumor_gt = (dgt/STAR_DIINDEL_GRID::SIZE);
}

struct is_nonsom_maker_t
{
    is_nonsom_maker_t();

    std::vector<bool> val;
};
}

std::ostream& operator<<(std::ostream& os,const DDIINDEL_GRID::index_t dgt);



// object used to pre-compute priors:
struct somatic_indel_caller_grid : private boost::noncopyable
{

    somatic_indel_caller_grid(const strelka_options& opt,
                              const indel_digt_caller& in_caller);

    void
    get_somatic_indel(const strelka_options& opt,
                      const strelka_deriv_options& dopt,
                      const starling_sample_options& normal_opt,
                      const starling_sample_options& tumor_opt,
                      const double indel_error_prob,
                      const double ref_error_prob,
                      const indel_key& ik,
                      const indel_data& normal_id,
                      const indel_data& tumor_id,
                      const bool is_use_alt_indel,
                      somatic_indel_call& sindel) const;

    void
    set_normal_prior(std::vector<blt_float_t>& normal_prior,
                     const double ref_error_prob,
                     const strelka_options& opt) const;

private:
    struct prior_set
    {
        prior_set()
            : normal(STAR_DIINDEL_GRID::SIZE)
            , normal_poly(STAR_DIINDEL_GRID::SIZE)
        {}

        std::vector<blt_float_t> normal;
        std::vector<blt_float_t> normal_poly;
    };

    prior_set _bare_lnprior;

    double _ln_som_match;
    double _ln_som_mismatch;
};
