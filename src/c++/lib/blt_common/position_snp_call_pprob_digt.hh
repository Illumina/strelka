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

#include "blt_common/blt_shared.hh"
#include "blt_common/snp_pos_info.hh"

#include "blt_util/digt.hh"
#include "blt_util/qscore.hh"

#include <boost/utility.hpp>

#include <array>
#include <iosfwd>



struct diploid_genotype
{
    diploid_genotype()
      : phredLoghood(DIGT::SIZE,0)
    {
        reset();
    }

    bool
    is_haploid() const
    {
        return (ploidy==1);
    }

    bool
    is_noploid() const
    {
        return (ploidy==0);
    }

    void reset()
    {
        is_snp=false;
        ploidy=2;
        ref_gt=0;
        genome.reset();
        poly.reset();
        std::fill(phredLoghood.begin(), phredLoghood.end(), 0);
    }

    struct result_set
    {
        result_set()
        {
            reset();
        }

        void
        reset()
        {
            max_gt=0;
            static const double p(1./static_cast<double>(pprob.size()));
            static const int qp(error_prob_to_qphred((1.-p)));
            snp_qphred=qp;
            max_gt_qphred=qp;
            std::fill(pprob.begin(),pprob.end(),p);
        }

        unsigned max_gt;
        int snp_qphred;
        int max_gt_qphred;
        std::array<double,DIGT::SIZE> pprob; // note this is intentionally stored at higher float resolution than the rest of the computation
    };

    // only used for PLs
    static const int maxQ;

    bool is_snp;

    /// a cheap way to add haploid calling capability, better solution: either haploid calls have their own object
    /// or this object is generalized to any ploidy
    ///
    int ploidy;

    unsigned ref_gt;
    result_set genome;
    result_set poly;

    std::vector<unsigned> phredLoghood;
};


// debuging output -- produces labels
//
std::ostream& operator<<(std::ostream& os,const diploid_genotype& dgt);


// more debuging output:
void
debug_dump_digt_lhood(const blt_float_t* lhood,
                      std::ostream& os);



// Use caller object to precalculate prior distributions based on
// theta value:
//
struct pprob_digt_caller : private boost::noncopyable
{
    explicit
    pprob_digt_caller(
        const blt_float_t theta);

    /// \brief call a snp @ pos by calculating the posterior probability
    /// of all possible genotypes for a diploid individual.
    ///
    /// When is_always_test is true, probabilities are calculated even
    /// when a snp could not exist at the site.
    ///
    void
    position_snp_call_pprob_digt(
        const blt_options& opt,
        const extended_pos_info& epi,
        diploid_genotype& dgt,
        double& strand_bias,
        const bool is_always_test = false) const;


    const blt_float_t*
    lnprior_genomic(
        const unsigned ref_id,
        const bool is_haploid = false) const
    {
        return get_prior(is_haploid)[ref_id].genome;
    }

    const blt_float_t*
    lnprior_polymorphic(
        const unsigned ref_id,
        const bool is_haploid = false) const
    {
        return get_prior(is_haploid)[ref_id].poly;
    }

    static
    void
    get_diploid_gt_lhood(
        const blt_options& opt,
        const extended_pos_info& epi,
        const bool is_het_bias,
        const blt_float_t het_bias,
        blt_float_t* const lhood,
        const bool is_strand_specific = false,
        const bool is_ss_fwd = false);

    static
    void
    calculate_result_set(
        const blt_float_t* lhood,
        const blt_float_t* lnprior,
        const unsigned ref_gt,
        diploid_genotype::result_set& rs);

    struct prior_set
    {
        prior_set()
        {
            std::fill(genome,genome+DIGT::SIZE,0);
            std::fill(poly,poly+DIGT::SIZE,0);
        }

        blt_float_t genome[DIGT::SIZE];
        blt_float_t poly[DIGT::SIZE];
    };

    typedef std::array<prior_set,N_BASE+1> prior_group;

private:

    const prior_group&
    get_prior(
        const bool is_haploid) const
    {
        return (is_haploid ? _lnprior_haploid : _lnprior);
    }

    prior_group _lnprior;
    prior_group _lnprior_haploid;
};
