#pragma once

#include "somatic_result_set.hh"
#include "strelka_shared.hh"

#include "blt_util/math_util.hh"
#include "blt_common/position_snp_call_pprob_digt.hh"

#include <cmath>

namespace SOMATIC_DIGT
{
enum index_t
{
    REF,
    HOM,
    HET,
    SIZE
};

inline
const char*
label(const unsigned idx)
{
    switch (idx)
    {
    case REF:
        return "ref";
    case HOM:
        return "hom";
    case HET:
        return "het";
    default:
        return "xxx";
    }
}

}

namespace TWO_STATE_SOMATIC
{
enum index_t
{
    NON_SOMATIC,
    SOMATIC,
    SIZE
};

}

namespace DIGT_GRID
{

// HET_RES is the number of points sampled between 0 and 0.5 on
// the continuous frequency scale. Thus a fully sampled axis will
// be sampled HET_RES*2+3 times.
//

//// The first three states are designed to overlap with
//// STAR_DIINDEL (ie. conventional diploid indel model), after this
//// the grid frequency (ie. approximations of continuous frequency)
//// states are added. The grid states are treated just like the
//// STAR_DIINDEL het state for certain purposes (printing, for
//// instance)

enum constants { HET_RES = 9,
                 HET_COUNT = HET_RES*2+1,
                 HOM_SIZE = 2,
                 PRESTRAND_SIZE = HOM_SIZE+HET_COUNT,
                 STRAND_STATE_SIZE = HET_RES
               };

enum index_t { SIZE = PRESTRAND_SIZE+STRAND_STATE_SIZE };

}

namespace DDIGT
{

enum index_t { SIZE = SOMATIC_DIGT::SIZE*TWO_STATE_SOMATIC::SIZE };

inline
unsigned
get_state(
    const unsigned normal_gt,
    const unsigned tumor_gt)
{
    return normal_gt*TWO_STATE_SOMATIC::SIZE + tumor_gt;
}

void
write_indel_state(const DDIGT::index_t dgt,
            std::ostream& os);

void
write_snv_state(const DDIGT::index_t dgt,
            const char ref_base,
            const char normal_alt_base,
            const char tumor_alt_base,
            std::ostream& os);

inline
void
get_digt_states(
    const unsigned dgt,
    unsigned& normal_gt,
    unsigned& tumor_gt)
{
    normal_gt = (dgt/TWO_STATE_SOMATIC::SIZE);
    tumor_gt = (dgt%TWO_STATE_SOMATIC::SIZE);
}

void
write_alt_alleles(unsigned alt_gt,
                  std::ostream& os);

}

namespace DDIGT_GRID
{

enum constants { PRESTRAND_SIZE = DIGT_GRID::PRESTRAND_SIZE*DIGT_GRID::PRESTRAND_SIZE };
enum index_t { SIZE = PRESTRAND_SIZE+DIGT_GRID::STRAND_STATE_SIZE };

inline
unsigned
get_state(
    const unsigned normal_gt,
    const unsigned tumor_gt)
{
    if (normal_gt<DIGT_GRID::PRESTRAND_SIZE) return normal_gt+DIGT_GRID::PRESTRAND_SIZE*tumor_gt;
    return PRESTRAND_SIZE+normal_gt-DIGT_GRID::PRESTRAND_SIZE;
}

}


void
calculate_bare_lnprior(const double theta,
        blt_float_t *bare_lnprior);

void
set_prior(
        const blt_float_t ssnv_freq_ratio,
        const blt_float_t ln_se_rate,   // ln (somatic error rate)
        const blt_float_t ln_cse_rate,  // ln (1 - somatic_error_rate)
        std::vector<blt_float_t>& ln_freq_given_somatic
        );

void
calculate_result_set_grid(
        const blt_float_t* normal_lhood,
        const blt_float_t* tumor_lhood,
        const std::vector<blt_float_t>& ln_freq_given_somatic,
        const blt_float_t* bare_lnprior_normal,
        const blt_float_t lnmatch,
        const blt_float_t lnmismatch,
        result_set& rs
        );
