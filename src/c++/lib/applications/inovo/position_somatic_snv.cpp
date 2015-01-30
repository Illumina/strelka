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
#include "position_somatic_snv.hh"

#include "blt_common/snp_util.hh"
#include "blt_util/log.hh"
#include "blt_util/prob_util.hh"
#include "blt_util/seq_util.hh"

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <iostream>
#include <iomanip>



std::ostream&
operator<<(std::ostream& os,
           const DDIGT::index_t dgt)
{

    unsigned normal_gt;
    unsigned tumor_gt;
    DDIGT::get_digt_states(dgt,normal_gt,tumor_gt);

    os << DIGT::label(normal_gt)
       << "->"
       << DIGT::label(tumor_gt);

    return os;

}



static
void
get_prior(const blt_float_t lnmatch,
          const blt_float_t lnmismatch,
          const blt_float_t* normal_lnprior,
          blt_float_t* prior)
{

    for (unsigned ngt(0); ngt<DIGT::SIZE; ++ngt)
    {
        const blt_float_t gtlnmatch(normal_lnprior[ngt]+lnmatch);
        const blt_float_t gtlnmismatch(normal_lnprior[ngt]+lnmismatch);
        for (unsigned tgt(0); tgt<DIGT::SIZE; ++tgt)
        {
            prior[DDIGT::get_state(ngt,tgt)] = ((tgt==ngt) ? gtlnmatch : gtlnmismatch);
        }
    }
}



somatic_snv_caller::
somatic_snv_caller(const inovo_options& opt,
                   const pprob_digt_caller& pd_caller) : _opt(opt)
{

    const blt_float_t lnmatch(std::log(1.-opt.somatic_snv_rate));
    const blt_float_t lnmismatch(std::log(opt.somatic_snv_rate/(static_cast<blt_float_t>(DIGT::SIZE-1))));

    for (unsigned i(0); i<(N_BASE+1); ++i)
    {
        prior_set& ps(_lnprior[i]);
        std::fill(ps.genome,ps.genome+DDIGT::SIZE,0);
        std::fill(ps.poly,ps.poly+DDIGT::SIZE,0);
    }

    for (unsigned i(0); i<(N_BASE+1); ++i)
    {
        prior_set& ps(_lnprior[i]);
        get_prior(lnmatch,lnmismatch,pd_caller.lnprior_genomic(i),ps.genome);
        get_prior(lnmatch,lnmismatch,pd_caller.lnprior_polymorphic(i),ps.poly);
    }
}


namespace DIGT
{

inline
bool
is_loh(const unsigned ngt,
       const unsigned tgt)
{

    if ((! is_het(ngt)) || is_het(tgt)) return false;
    return (expect2(tgt,ngt)>0);
}
}



typedef somatic_snv_genotype::result_set result_set;



#ifdef SOMATIC_DEBUG
static
void
debug_dump_ddigt_lhood(const blt_float_t* lhood,
                       std::ostream& os)
{

    blt_float_t pprob[DDIGT::SIZE];
    for (unsigned gt(0); gt<DDIGT::SIZE; ++gt)
    {
        pprob[gt] = lhood[gt];
    }

    unsigned max_gt(0);
    normalize_ln_prob(pprob,pprob+DDIGT::SIZE,max_gt);

    os << std::setprecision(3) << std::fixed;
    for (unsigned ngt(0); ngt<DIGT::SIZE; ++ngt)
    {
        for (unsigned tgt(0); tgt<DIGT::SIZE; ++tgt)
        {
            const unsigned dgt(DDIGT::get_state(ngt,tgt));
            os << static_cast<DDIGT::index_t>(dgt) << ": " << -std::log(pprob[dgt]) << " ";
        }
        os << "\n";
    }
    os.unsetf(std::ios::fixed);
}
#endif



static
void
calculate_result_set(const blt_float_t* normal_lhood,
                     const blt_float_t* tumor_lhood,
                     const blt_float_t* lnprior,
                     const unsigned ref_gt,
                     result_set& rs)
{

    // mult by prior distro to get unnormalized pprob:
    //
    for (unsigned ngt(0); ngt<DIGT::SIZE; ++ngt)
    {
        for (unsigned tgt(0); tgt<DIGT::SIZE; ++tgt)
        {
            const unsigned dgt(DDIGT::get_state(ngt,tgt));
            rs.pprob[dgt] = normal_lhood[ngt]+tumor_lhood[tgt]+lnprior[dgt];
        }
    }

    normalize_ln_distro(rs.pprob.begin(),rs.pprob.end(),rs.max_gt);

    blt_float_t nonsomatic_sum(0);
    for (unsigned gt(0); gt<DIGT::SIZE; ++gt)
    {
        nonsomatic_sum += rs.pprob[DDIGT::get_state(gt,gt)];
    }
    rs.snv_qphred=error_prob_to_qphred(nonsomatic_sum);

    blt_float_t not_somfrom_ref_sum(nonsomatic_sum);
    blt_float_t not_somfrom_het_sum(nonsomatic_sum);
    blt_float_t not_somfrom_het_loh_sum(nonsomatic_sum);
    blt_float_t not_somfrom_het_nonloh_sum(nonsomatic_sum);
    blt_float_t not_somfrom_hom_sum(nonsomatic_sum);
    blt_float_t not_somfrom_anyhom_sum(nonsomatic_sum);
    for (unsigned ngt(0); ngt<DIGT::SIZE; ++ngt)
    {
        const bool is_ref(ref_gt == ngt);
        const bool is_het(DIGT::is_het(ngt));
        for (unsigned tgt(0); tgt<DIGT::SIZE; ++tgt)
        {
            if (ngt==tgt) continue;
            const blt_float_t val(rs.pprob[DDIGT::get_state(ngt,tgt)]);
            if (not is_ref) not_somfrom_ref_sum += val;
            if (is_het)
            {
                not_somfrom_anyhom_sum += val;

                const bool is_loh(DIGT::is_loh(ngt,tgt));
                if (is_loh) not_somfrom_het_nonloh_sum += val;
                else       not_somfrom_het_loh_sum += val;
            }
            else
            {
                not_somfrom_het_sum += val;
                not_somfrom_het_loh_sum += val;
                not_somfrom_het_nonloh_sum += val;
            }
            if (is_het or is_ref) not_somfrom_hom_sum += val;
        }
    }
    rs.snv_from_ref_qphred=error_prob_to_qphred(not_somfrom_ref_sum);
    rs.snv_from_het_qphred=error_prob_to_qphred(not_somfrom_het_sum);
    rs.snv_from_het_loh_qphred=error_prob_to_qphred(not_somfrom_het_loh_sum);
    rs.snv_from_het_nonloh_qphred=error_prob_to_qphred(not_somfrom_het_nonloh_sum);
    rs.snv_from_hom_qphred=error_prob_to_qphred(not_somfrom_hom_sum);
    rs.snv_from_anyhom_qphred=error_prob_to_qphred(not_somfrom_anyhom_sum);

    rs.max_gt_qphred=error_prob_to_qphred(prob_comp(rs.pprob.begin(),rs.pprob.end(),rs.max_gt));
}



///
///
void
somatic_snv_caller::
position_somatic_snv_call(const extended_pos_info& normal_epi,
                          const extended_pos_info& tumor_epi,
                          const extended_pos_info* normal_epi_t2_ptr,
                          const extended_pos_info* tumor_epi_t2_ptr,
                          somatic_snv_genotype& sgt) const
{

    static const bool is_always_test(false);

    {
        const snp_pos_info& normal_pi(normal_epi.pi);
        const snp_pos_info& tumor_pi(tumor_epi.pi);

        if (normal_pi.ref_base=='N') return;
        sgt.ref_gt=base_to_id(normal_pi.ref_base);

        // check that a non-reference call meeting quality criteria even
        // exists:
        if (! is_always_test)
        {
            if (is_spi_allref(normal_pi,sgt.ref_gt) && is_spi_allref(tumor_pi,sgt.ref_gt)) return;
        }
    }

    // strawman model treats normal and tumor as independent, so
    // calculate separate lhoods:
    blt_float_t normal_lhood[DIGT::SIZE];
    blt_float_t tumor_lhood[DIGT::SIZE];

    const bool is_tier2(NULL != normal_epi_t2_ptr);

    static const unsigned n_tier(2);
    result_set tier_rs[n_tier];
    for (unsigned i(0); i<n_tier; ++i)
    {
        const bool is_include_tier2(i==1);
        if (is_include_tier2)
        {
            if (not is_tier2) continue;
            if (tier_rs[0].snv_qphred==0)
            {
                tier_rs[1].snv_qphred=0;
                continue;
            }
        }

        // get likelihood of each genotype
        //
        static const bool is_normal_het_bias(false);
        static const blt_float_t normal_het_bias(0.0);
        static const bool is_tumor_het_bias(false);
        static const blt_float_t tumor_het_bias(0.0);

        const extended_pos_info& nepi(is_include_tier2 ? *normal_epi_t2_ptr : normal_epi );
        const extended_pos_info& tepi(is_include_tier2 ? *tumor_epi_t2_ptr : tumor_epi );
        pprob_digt_caller::get_diploid_gt_lhood(_opt,nepi,is_normal_het_bias,normal_het_bias,normal_lhood);
        pprob_digt_caller::get_diploid_gt_lhood(_opt,tepi,is_tumor_het_bias,tumor_het_bias,tumor_lhood);

        // genomic site results:
        calculate_result_set(normal_lhood,tumor_lhood,
                             lnprior_genomic(sgt.ref_gt),sgt.ref_gt,tier_rs[i]);

#ifdef ENABLE_POLY
        // polymorphic site results:
        assert(0); // still needs to be adapted for 2-tier system:
        calculate_result_set(normal_lhood,tumor_lhood,
                             lnprior_polymorphic(sgt.ref_gt),sgt.ref_gt,sgt.poly);
#else
        sgt.poly.snv_qphred = 0;
#endif
    }

    sgt.tier=0;
    if (is_tier2)
    {
        if (tier_rs[0].snv_qphred > tier_rs[1].snv_qphred)
        {
            sgt.tier=1;
        }
    }
    sgt.genome=tier_rs[sgt.tier];


    // not sure which prior to use yet:
    //
    sgt.is_snv=((sgt.genome.snv_qphred != 0) or (sgt.poly.snv_qphred != 0));


#ifdef SOMATIC_DEBUG
    if (sgt.is_snv)
    {
        log_os << "tier1_qphred: " << tier_rs[0].snv_qphred << "\n";
        log_os << "tier2_qphred: " << tier_rs[1].snv_qphred << "\n";
        result_set rs(sgt.genome);
        const blt_float_t* lnprior(lnprior_genomic(sgt.ref_gt));
        for (unsigned ngt(0); ngt<DIGT::SIZE; ++ngt)
        {
            for (unsigned tgt(0); tgt<DIGT::SIZE; ++tgt)
            {
                const unsigned dgt(DDIGT::get_state(ngt,tgt));
                rs.pprob[dgt] = normal_lhood[ngt]+tumor_lhood[tgt]+lnprior[dgt];
            }
        }

        log_os << "normal_digt: ";
        debug_dump_digt_lhood(normal_lhood,log_os);
        log_os << "\n";
        log_os << "tumor_digt: ";
        debug_dump_digt_lhood(tumor_lhood,log_os);
        log_os << "\n";
        log_os << "prior:\n";
        debug_dump_ddigt_lhood(lnprior,log_os);
        log_os << "pprob:\n";
        blt_float_t pprob[DDIGT::SIZE];
        for (unsigned gt(0); gt<DDIGT::SIZE; ++gt)
        {
            pprob[gt] = rs.pprob[gt];
        }
        debug_dump_ddigt_lhood(pprob,log_os);
    }
#endif
}



static
void
write_result_set(const result_set& rs,
                 std::ostream& os)
{
    os << rs.snv_qphred
       << '\t' << rs.snv_from_ref_qphred
       << '\t' << rs.snv_from_het_qphred
       << '\t' << rs.snv_from_het_loh_qphred
       << '\t' << rs.snv_from_het_nonloh_qphred
       << '\t' << rs.snv_from_hom_qphred
       << '\t' << rs.snv_from_anyhom_qphred
       << '\t' << static_cast<DDIGT::index_t>(rs.max_gt)
       << '\t' << rs.max_gt_qphred;
}



void
write_somatic_snv_genotype(const inovo_options& opt,
                           const somatic_snv_genotype& sgt,
                           const snp_pos_info& normal_pi,
                           const snp_pos_info& tumor_pi,
                           std::ostream& os)
{
    os << std::setprecision(10) << std::fixed;

    const result_set& ge(sgt.genome);
#ifdef ENABLE_POLY
    const result_set& po(sgt.poly);
#endif

    os << (sgt.tier+1) << '\t';
    write_result_set(ge,os);
#ifdef ENABLE_POLY
    os << '\t';
    write_result_set(po,os);
#endif
    if (opt.is_print_used_allele_counts)
    {
        normal_pi.print_known_counts(os,opt.used_allele_count_min_qscore);
        tumor_pi.print_known_counts(os,opt.used_allele_count_min_qscore);
    }

    os.unsetf(std::ios::fixed);
}
