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

#if 0

#include "somatic_indel.hh"

#include "blt_util/math_util.hh"
#include "blt_util/prob_util.hh"

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <iostream>
#include <iomanip>


#ifdef DEBUG_INDEL_CALL
#include "blt_util/log.hh"
#endif



std::ostream&
operator<<(std::ostream& os,
           const DDIINDEL::index_t dgt) {

    unsigned normal_gt;
    unsigned tumor_gt;
    DDIINDEL::get_digt_states(dgt,normal_gt,tumor_gt);

    os << STAR_DIINDEL::label(normal_gt)
       << "->"
       << STAR_DIINDEL::label(tumor_gt);

    return os;

}



somatic_indel_caller::
somatic_indel_caller(const strelka_options& opt,
                     const indel_digt_caller& in_caller)
    : _in_caller(in_caller)
{
    _ln_som_match=(log1p_switch(-opt.somatic_indel_rate));
    _ln_som_mismatch=(std::log(opt.somatic_indel_rate/(static_cast<double>(STAR_DIINDEL::SIZE-1))));
}



typedef somatic_indel_call::result_set result_set;



static
void
calculate_result_set(const strelka_options& opt,
                     const double* normal_lnprior,
                     const double lnmatch,
                     const double lnmismatch,
                     const double* normal_lhood,
                     const double* tumor_lhood,
                     result_set& rs) {

#ifdef SOMATIC_DEBUG
    std::vector<double> check_prior(DDIINDEL::SIZE);

    for(unsigned ngt(0); ngt<STAR_DIINDEL::SIZE; ++ngt) {
        const double base_prior(normal_lnprior[ngt]);
        for(unsigned tgt(0); tgt<STAR_DIINDEL::SIZE; ++tgt) {
            const unsigned dgt(DDIINDEL::get_state(ngt,tgt));
            check_prior[dgt] =
                base_prior+
                ((tgt==ngt) ? lnmatch : lnmismatch);
        }
    }

    check_ln_distro(check_prior.begin(),
                    check_prior.end(),
                    "somatic indel full prior");
#endif

    // get unnormalized posterior:
    std::vector<double> pprob(DDIINDEL::SIZE);

    for(unsigned ngt(0); ngt<STAR_DIINDEL::SIZE; ++ngt) {
        const double base_prior(normal_lnprior[ngt]);
        for(unsigned tgt(0); tgt<STAR_DIINDEL::SIZE; ++tgt) {
            const unsigned dgt(DDIINDEL::get_state(ngt,tgt));
            pprob[dgt] =
                normal_lhood[ngt]+
                tumor_lhood[tgt]+
                base_prior+
                ((tgt==ngt) ? lnmatch : lnmismatch);
        }
    }

    normalize_ln_distro(pprob.begin(),pprob.end(),rs.max_gt);

#ifdef DEBUG_INDEL_CALL
    log_os << "INDEL_CALL pprob(noindel),pprob(hom),pprob(het): " << pprob[STAR_DIINDEL::NOINDEL] << " " << pprob[STAR_DIINDEL::HOM] << " " << pprob[STAR_DIINDEL::HET] << "\n";
#endif
    double nonsomatic_sum(0);
    for(unsigned gt(0); gt<STAR_DIINDEL::SIZE; ++gt) {
        nonsomatic_sum += pprob[DDIINDEL::get_state(gt,gt)];
    }
    rs.sindel_qphred=error_prob_to_qphred(nonsomatic_sum);

    double not_somfrom_sum[STAR_DIINDEL::SIZE];
    for(unsigned sgt(0); sgt<STAR_DIINDEL::SIZE; ++sgt) {
        not_somfrom_sum[sgt]=nonsomatic_sum;
        for(unsigned ngt(0); ngt<STAR_DIINDEL::SIZE; ++ngt) {
            if(sgt==ngt) continue;
            for(unsigned tgt(0); tgt<STAR_DIINDEL::SIZE; ++tgt) {
                if(tgt==ngt) continue;
                not_somfrom_sum[sgt] += pprob[DDIINDEL::get_state(ngt,tgt)];
            }
        }
    }
    rs.sindel_from_ref_qphred=error_prob_to_qphred(not_somfrom_sum[STAR_DIINDEL::NOINDEL]);
    rs.sindel_from_het_qphred=error_prob_to_qphred(not_somfrom_sum[STAR_DIINDEL::HET]);
    rs.sindel_from_hom_qphred=error_prob_to_qphred(not_somfrom_sum[STAR_DIINDEL::HOM]);

    double not_somfromanyhom_sum(nonsomatic_sum);
    for(unsigned ngt(0); ngt<STAR_DIINDEL::SIZE; ++ngt) {
        if(STAR_DIINDEL::HET != ngt) continue;
        for(unsigned tgt(0); tgt<STAR_DIINDEL::SIZE; ++tgt) {
            if(tgt==ngt) continue;
            not_somfromanyhom_sum += pprob[DDIINDEL::get_state(ngt,tgt)];
        }
    }
    rs.sindel_from_anyhom_qphred=error_prob_to_qphred(not_somfromanyhom_sum);

    rs.max_gt_qphred=error_prob_to_qphred(prob_comp(pprob.begin(),pprob.end(),rs.max_gt));
}



///
///
void
somatic_indel_caller::
get_somatic_indel(const strelka_options& opt,
                  const strelka_deriv_options& dopt,
                  const double indel_error_prob,
                  const double ref_error_prob,
                  const indel_key& ik,
                  const indel_data& normal_id,
                  const indel_data& tumor_id,
                  const bool is_use_alt_indel,
                  somatic_indel_call& sindel) const {

    // for now, lhood calculation of tumor and normal are independent:

    // get likelihood of each genotype
    static const bool is_normal_het_bias(false);
    static const bool is_tumor_het_bias(true);
    static const double normal_het_bias(0.0);
    static const double tumor_het_bias(0.4);
    double normal_lhood[STAR_DIINDEL::SIZE];
    double tumor_lhood[STAR_DIINDEL::SIZE];

    static const unsigned n_tier(2);
    result_set tier_rs[n_tier];
    for(unsigned i(0); i<n_tier; ++i) {
        const bool is_include_tier2(i==1);
        if(is_include_tier2) {
            if(not opt.is_tier2()) continue;
            if(tier_rs[0].sindel_qphred==0) {
                tier_rs[1].sindel_qphred=0;
                continue;
            }
        }
        indel_digt_caller::get_indel_digt_lhood(opt,dopt,indel_error_prob,ref_error_prob,ik,normal_id,
                                                is_normal_het_bias,normal_het_bias,
                                                is_include_tier2,is_use_alt_indel,
                                                normal_lhood);
        indel_digt_caller::get_indel_digt_lhood(opt,dopt,indel_error_prob,ref_error_prob,ik,tumor_id,
                                                is_tumor_het_bias,tumor_het_bias,
                                                is_include_tier2,is_use_alt_indel,
                                                tumor_lhood);

        calculate_result_set(opt,normal_lnprior(),
                             _ln_som_match,_ln_som_mismatch,
                             normal_lhood,tumor_lhood,tier_rs[i]);
    }

    sindel.tier=0;
    if(opt.is_tier2()) {
        if(tier_rs[0].sindel_qphred > tier_rs[1].sindel_qphred) {
            sindel.tier=1;
        }
    }

    sindel.rs=tier_rs[sindel.tier];
    sindel.is_indel=(sindel.rs.sindel_qphred != 0);
}



static
void
write_isri(const starling_indel_sample_report_info& isri,
           std::ostream& os) {

    os << '\t' << isri.depth
       << '\t' << isri.n_q30_ref_reads+isri.n_q30_alt_reads
       << '\t' << isri.n_q30_indel_reads
       << '\t' << isri.n_other_reads;
}



void
write_somatic_indel_file(const somatic_indel_call& sindel,
                         const starling_indel_report_info& iri,
                         const starling_indel_sample_report_info& nisri,
                         const starling_indel_sample_report_info& tisri,
                         std::ostream& os) {

    const somatic_indel_call::result_set& rs(sindel.rs);

    os << iri.desc
       << '\t' << iri.ref_upstream
       << '\t' << iri.ref_seq << "/" << iri.indel_seq
       << '\t' << iri.ref_downstream
       << '\t' << (sindel.tier+1)
       << '\t' << rs.sindel_qphred
       << '\t' << rs.sindel_from_ref_qphred
       << '\t' << rs.sindel_from_het_qphred
       << '\t' << rs.sindel_from_hom_qphred
       << '\t' << rs.sindel_from_anyhom_qphred
       << '\t' << static_cast<DDIINDEL::index_t>(rs.max_gt)
       << '\t' << rs.max_gt_qphred;

    write_isri(nisri,os);
    write_isri(tisri,os);

    os << '\t' << iri.repeat_unit
       << '\t' << iri.ref_repeat_count
       << '\t' << iri.indel_repeat_count;
}


#endif
