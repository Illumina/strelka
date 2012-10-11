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

/// \author Chris Saunders
///
#ifndef __POSITION_SNP_CALL_PPROB_DIGT_HH
#define __POSITION_SNP_CALL_PPROB_DIGT_HH

#include "blt_common/blt_shared.hh"
#include "blt_common/snp_pos_info.hh"

#include "blt_util/digt.hh"
#include "blt_util/qscore.hh"

#include <boost/array.hpp>
#include <boost/utility.hpp>

#include <iosfwd>



struct diploid_genotype {

    diploid_genotype() 
        : is_snp(false), ref_gt(0) {}

    struct result_set {

        result_set() 
            : max_gt(0)
        {
            static const blt_float_t p(1./static_cast<blt_float_t>(DIGT::SIZE));
            static const int qp(error_prob_to_qphred((1.-p)));
            snp_qphred=qp;
            max_gt_qphred=qp;
            for(unsigned i(0);i<DIGT::SIZE;++i) {
                pprob[i] = p;
            } 
        }
  
        unsigned max_gt;
        int snp_qphred;
        int max_gt_qphred;
        boost::array<double,DIGT::SIZE> pprob; // note this is intentionally stored at higher float resolution than the rest of the computation
    };

    bool is_snp;
    unsigned ref_gt;
    result_set genome;
    result_set poly;
    double sb;
};


// debuging output -- produces labels
//
std::ostream& operator<<(std::ostream& os,const diploid_genotype& dgt);


// more debuging output:
void
debug_dump_digt_lhood(const blt_float_t* lhood,
                      std::ostream& os);


// allele call output:
//
void
write_diploid_genotype_allele(const blt_options& opt,
                              const snp_pos_info& pi,
                              const diploid_genotype& dgt,
                              std::ostream& os,
                              const unsigned hpol);

// snp call output:
//
inline
void
write_diploid_genotype_snp(const blt_options& opt,
                           const snp_pos_info& pi,
                           const diploid_genotype& dgt,
                           std::ostream& os,
                           const unsigned hpol){

    write_diploid_genotype_allele(opt,pi,dgt,os,hpol);
}


// Use caller object to precalculate prior distributions based on
// theta value:
//
struct pprob_digt_caller : private boost::noncopyable {

    pprob_digt_caller(const blt_float_t theta);

    /// \brief call a snp @ pos by calculating the posterior probability
    /// of all possible genotypes for a diploid individual.
    ///
    /// When is_always_test is true, probabilities are calculated even
    /// when a snp could not exist at the site.
    ///
    void
    position_snp_call_pprob_digt(const blt_options& opt,
                                 const extended_pos_info& epi,
                                 diploid_genotype& dgt,
                                 const bool is_always_test = false) const;


    const blt_float_t*
    lnprior_genomic(const unsigned ref_id) const {
        return _lnprior[ref_id].genome;
    }

    const blt_float_t*
    lnprior_polymorphic(const unsigned ref_id) const {
        return _lnprior[ref_id].poly;
    }

    static
    void
    get_diploid_gt_lhood(const blt_options& opt,
                         const extended_pos_info& epi,
                         const bool is_het_bias,
                         const blt_float_t het_bias,
                         blt_float_t* const lhood,
                         const bool is_strand_specific = false,
                         const bool is_ss_fwd = false);

    static
    void
    calculate_result_set(const blt_float_t* lhood,
                         const blt_float_t* lnprior,
                         const unsigned ref_gt,
                         diploid_genotype::result_set& rs);

private:
    struct prior_set {
        blt_float_t genome[DIGT::SIZE];
        blt_float_t poly[DIGT::SIZE];
    };

    prior_set _lnprior[N_BASE+1];
};

#endif
