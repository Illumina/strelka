//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
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
        ploidy=2;
        ref_gt=0;
        strand_bias=0;
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
            ref_pprob=0.;
            snp_qphred=0;
            max_gt_qphred=0;
        }

        unsigned max_gt;
        double ref_pprob;
        int snp_qphred;
        int max_gt_qphred;
    };

    // only used for PLs
    static const int maxQ;

    bool is_snp() const
    {
        return (genome.snp_qphred != 0);
    }

    /// a cheap way to add haploid calling capability, better solution: either haploid calls have their own object
    /// or this object is generalized to any ploidy
    ///
    int ploidy;

    unsigned ref_gt;
    double strand_bias;
    result_set genome;
    result_set poly;
    std::vector<unsigned> phredLoghood;
};


// debugging output -- produces labels
//
std::ostream& operator<<(std::ostream& os,const diploid_genotype& dgt);


// more debugging output:
void
debug_dump_digt_lhood(const blt_float_t* lhood,
                      std::ostream& os);



/// Used to precalculate prior distributions based on theta value
///
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
        const bool useHetVariantFrequencyExtension,
        const blt_float_t hetVariantFrequencyExtension,
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
