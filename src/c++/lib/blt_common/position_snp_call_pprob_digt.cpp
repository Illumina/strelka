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

#include "blt_common/position_snp_call_pprob_digt.hh"
#include "blt_common/snp_util.hh"
#include "blt_util/log.hh"
#include "blt_util/logSumUtil.hh"
#include "blt_util/prob_util.hh"
#include "blt_util/seq_util.hh"

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <iostream>
#include <iomanip>



const int diploid_genotype::maxQ = 999;


const blt_float_t one_third(1./3.);
const blt_float_t log_one_third(std::log(one_third));
const blt_float_t one_half(1./2.);
const blt_float_t log_one_half(std::log(one_half));



static
void
get_genomic_prior(
    const unsigned ref_gt,
    const blt_float_t theta,
    blt_float_t* const prior)
{
    blt_float_t prior_sum(0.);
    for (unsigned gt(0); gt<DIGT::SIZE; ++gt)
    {
        if (gt==ref_gt) continue;
        prior[gt]=(theta*one_third);
        if (DIGT::is_het(gt))
        {
            if (DIGT::expect(ref_gt,gt)<=0.) prior[gt]*=theta;
        }
        else
        {
            prior[gt]*=.5;
        }
        prior_sum += prior[gt];
    }
    assert(prior_sum <= 1.);
    prior[ref_gt] = (1.-prior_sum);
}



static
void
get_haploid_genomic_prior(
    const unsigned ref_gt,
    const blt_float_t theta,
    blt_float_t* const prior)
{
    blt_float_t prior_sum(0.);
    for (unsigned gt(0); gt<DIGT::SIZE; ++gt)
    {
        if (gt==ref_gt) continue;
        if (DIGT::is_het(gt))
        {
            prior[gt]=0;
        }
        else
        {
            prior[gt]=(theta*one_third);
        }
        prior_sum += prior[gt];
    }
    assert(prior_sum <= 1.);
    prior[ref_gt] = (1.-prior_sum);
}



static
void
get_poly_prior(
    const unsigned ref_gt,
    const blt_float_t theta,
    blt_float_t* const prior)
{
    blt_float_t prior_sum(0.);
    const blt_float_t ctheta(1.-theta);
    for (unsigned gt(0); gt<DIGT::SIZE; ++gt)
    {
        if (gt==ref_gt)
        {
            prior[gt]=0.25*(ctheta);
        }
        else if (DIGT::is_het(gt))
        {
            if (DIGT::expect(ref_gt,gt)<=0.)
            {
                prior[gt] = theta*one_third;
            }
            else
            {
                prior[gt] = 0.5*one_third*ctheta;
            }
        }
        else
        {
            prior[gt] = 0.25*one_third*ctheta;
        }
        prior_sum += prior[gt];
    }
    assert(std::abs(1.-prior_sum) < 0.0001);
}



static
void
get_haploid_poly_prior(
    const unsigned ref_gt,
    const blt_float_t /*theta*/,
    blt_float_t* const prior)
{
    blt_float_t prior_sum(0.);
    for (unsigned gt(0); gt<DIGT::SIZE; ++gt)
    {
        if (gt==ref_gt)
        {
            prior[gt]=0.5;
        }
        else if (DIGT::is_het(gt))
        {
            prior[gt]=0;
        }
        else
        {
            prior[gt] = 0.5*one_third;
        }
        prior_sum += prior[gt];
    }
    assert(std::abs(1.-prior_sum) < 0.0001);
}



static
void
sum_gt(blt_float_t* const x1,
       const blt_float_t* const x2)
{
    for (unsigned gt(0); gt<DIGT::SIZE; ++gt)
    {
        x1[gt] += x2[gt];
    }
}



static
void
norm_gt(blt_float_t* const x)
{
    blt_float_t sum(0);
    for (unsigned gt(0); gt<DIGT::SIZE; ++gt)
    {
        sum += x[gt];
    }
    sum = 1./sum;
    for (unsigned gt(0); gt<DIGT::SIZE; ++gt)
    {
        x[gt] *= sum;
    }
}



static
void
log_gt(blt_float_t* const x)
{
    for (unsigned gt(0); gt<DIGT::SIZE; ++gt)
    {
        x[gt] = std::log(x[gt]);
    }
}



static
void
finish_prior(
    pprob_digt_caller::prior_group& prior)
{
    // 'N' prior is the average:
    auto& nps(prior[N_BASE]);
    for (unsigned i(0); i<N_BASE; ++i)
    {
        auto& ps(prior[i]);
        sum_gt(nps.genome,ps.genome);
        sum_gt(nps.poly,ps.poly);
    }
    norm_gt(nps.genome);
    norm_gt(nps.poly);

    // take logs:
    for (unsigned i(0); i<(N_BASE+1); ++i)
    {
        auto& ps(prior[i]);
        log_gt(ps.genome);
        log_gt(ps.poly);
    }
}



pprob_digt_caller::
pprob_digt_caller(
    const blt_float_t theta)
{
    for (unsigned i(0); i<N_BASE; ++i)
    {
        prior_set& ps(_lnprior[i]);
        get_genomic_prior(i,theta,ps.genome);
        get_poly_prior(i,theta,ps.poly);

        prior_set& psh(_lnprior_haploid[i]);
        get_haploid_genomic_prior(i,theta,psh.genome);
        get_haploid_poly_prior(i,theta,psh.poly);
    }

    finish_prior(_lnprior);
    finish_prior(_lnprior_haploid);
}



static
void
increment_het_ratio_lhood(const extended_pos_info& epi,
                          const blt_float_t het_ratio,
                          blt_float_t* all_het_lhood,
                          const bool is_strand_specific,
                          const bool is_ss_fwd)
{
    const blt_float_t chet_ratio(1-het_ratio);

    // multiply probs of alternate ratios into local likelihoods, then
    // *add* them to the global tally (effictively this is the sum lhood of
    // many different heterozygous genotypes).
    //
    // in the gt_high genotype, the first allele (in lexicographical
    // order) is expected at het_ratio and the second allele is
    // expected at chet_ratio.  gt_low genotype is vice versa.
    //
    blt_float_t lhood_high[DIGT::SIZE];
    blt_float_t lhood_low[DIGT::SIZE];
    for (unsigned gt(0); gt<DIGT::SIZE; ++gt)
    {
        lhood_high[gt] = 0.;
        lhood_low[gt] = 0.;
    }

    const snp_pos_info& pi(epi.pi);
    const unsigned ref_gt(base_to_id(pi.get_ref_base()));

    const unsigned n_calls(pi.calls.size());

    blt_float_t val_high[3];

    for (unsigned i(0); i<n_calls; ++i)
    {
        const blt_float_t eprob(epi.de[i]);
        const blt_float_t ceprob(1.-pi.calls[i].error_prob());

        // precalculate the result for expect values of 0.0, het_ratio, chet_ratio, 1.0
        val_high[0] = std::log(eprob)+log_one_third;
        val_high[1] = std::log((ceprob)*het_ratio+((1.-ceprob)*one_third)*chet_ratio);
        val_high[2] = std::log((ceprob)*chet_ratio+((1.-ceprob)*one_third)*het_ratio);

        const bool is_force_ref(is_strand_specific && (is_ss_fwd!=pi.calls[i].is_fwd_strand));

        const uint8_t obs_id(pi.calls[i].base_id);
        for (unsigned gt(N_BASE); gt<DIGT::SIZE; ++gt)
        {
            static const uint8_t low_remap[] = {0,2,1};
            const unsigned key(DIGT::expect2_bias(obs_id,(is_force_ref ? ref_gt : gt)));
            lhood_high[gt] += val_high[key];
            lhood_low[gt] += val_high[low_remap[key]];
        }
    }

    for (unsigned gt(0); gt<DIGT::SIZE; ++gt)
    {
        if (! DIGT::is_het(gt)) continue;
        all_het_lhood[gt] = getLogSum(all_het_lhood[gt], lhood_high[gt], lhood_low[gt]);
    }
}



void
pprob_digt_caller::
get_diploid_gt_lhood(const blt_options& opt,
                     const extended_pos_info& epi,
                     const bool useHetVariantFrequencyExtension,
                     const blt_float_t hetVariantFrequencyExtension,
                     blt_float_t* const lhood,
                     const bool is_strand_specific,
                     const bool is_ss_fwd)
{
    // get likelihood of each genotype
    for (unsigned gt(0); gt<DIGT::SIZE; ++gt) lhood[gt] = 0.;

    const snp_pos_info& pi(epi.pi);
    const unsigned ref_gt(base_to_id(pi.get_ref_base()));

    const unsigned n_calls(pi.calls.size());
    for (unsigned i(0); i<n_calls; ++i)
    {
        const base_call& bc(pi.calls[i]);
        const blt_float_t eprob(epi.de[i]);
        const blt_float_t ceprob(1.-bc.error_prob());
        const blt_float_t lnce(bc.ln_comp_error_prob());

        // precalculate the result for expect values of 0.0, 0.5 & 1.0
        blt_float_t val[3];
        val[0] = std::log(eprob)+log_one_third;
        val[1] = std::log((ceprob)+((1.-ceprob)*one_third))+log_one_half;
        val[2] = lnce;

        const bool is_force_ref(is_strand_specific && (is_ss_fwd!=bc.is_fwd_strand));

        const uint8_t obs_id(bc.base_id);
        for (unsigned gt(0); gt<DIGT::SIZE; ++gt)
        {
            lhood[gt] += val[DIGT::expect2(obs_id,(is_force_ref ? ref_gt : gt))];
        }
    }

    if (useHetVariantFrequencyExtension)
    {
        // loop is currently setup to assume a uniform het ratio subgenotype prior
        const unsigned n_bias_steps(1+static_cast<unsigned>(hetVariantFrequencyExtension/opt.maxHetVariantFrequencyIncrement));
        const blt_float_t ratio_increment(hetVariantFrequencyExtension/static_cast<blt_float_t>(n_bias_steps));
        for (unsigned i(0); i<n_bias_steps; ++i)
        {
            const blt_float_t het_ratio(0.5+(i+1)*ratio_increment);
            increment_het_ratio_lhood(epi,het_ratio,lhood,is_strand_specific,is_ss_fwd);
        }

        const unsigned n_het_subgt(1+2*n_bias_steps);
        const blt_float_t subgt_log_prior(-std::log(static_cast<blt_float_t>(n_het_subgt)));

        for (unsigned gt(0); gt<DIGT::SIZE; ++gt)
        {
            if (! DIGT::is_het(gt)) continue;
            lhood[gt] += subgt_log_prior;
        }
    }
}

typedef diploid_genotype::result_set result_set;

void
debug_dump_digt_lhood(const blt_float_t* lhood,
                      std::ostream& os)
{
    blt_float_t pprob[DIGT::SIZE];
    for (unsigned gt(0); gt<DIGT::SIZE; ++gt)
    {
        pprob[gt] = lhood[gt];
    }

    unsigned max_gt(0);
    normalizeLogDistro(pprob, pprob + DIGT::SIZE, max_gt);

    os << std::setprecision(3) << std::fixed;
    for (unsigned gt(0); gt<DIGT::SIZE; ++gt)
    {
        os << DIGT::label(gt) << ": " << -std::log(pprob[gt]) << " ";
    }
    os.unsetf(std::ios::fixed);
}



void
pprob_digt_caller::
calculate_result_set(const blt_float_t* lhood,
                     const blt_float_t* lnprior,
                     const unsigned ref_gt,
                     result_set& rs)
{
    std::array<double,DIGT::SIZE> pprob; // note this is intentionally stored at higher float resolution than the rest of the computation

    // mult by prior distro to get unnormalized pprob:
    //
    for (unsigned gt(0); gt<DIGT::SIZE; ++gt)
    {
        pprob[gt] = lhood[gt] + lnprior[gt];
    }

    normalizeLogDistro(pprob.begin(), pprob.end(), rs.max_gt);

    rs.ref_pprob=pprob[ref_gt];
    rs.snp_qphred=error_prob_to_qphred(pprob[ref_gt]);
    rs.max_gt_qphred=error_prob_to_qphred(prob_comp(pprob.begin(),pprob.end(),rs.max_gt));
}



///
/// \param is_always_test - continue the full computation even for a site which is
///                          obviously non-variant
///
///
/// The likelihood of an observed 'column' of observations for any
/// genotype: P(obs=ACT|genotype=AT) is a sum over all possible
/// columns:
///
/// sum { c in {A,C,G,T}**3 } { P(obs=ACT|col=c) * P(col=c|genotype=AT) }
///
/// Enumerating all {A,C,G,T}**col_size possible column values is
/// inefficent, so we only consider those columns for which
/// P(col|genotype) is nonzero.  For homozygous genotypes this always
/// leaves a single column; for heterozygous genotypes this is a sum
/// over all valid columns:
///
/// P(obs=ACT|genotype=AT)=
/// sum { c in {A,T}**3 } { P(obs=ACT|col=c) * P(col=c|genotype=AT }
///
/// that is...
///
/// P(obs|AAA)*P(AAA|genotype)+
/// P(obs|AAT)*P(AAT|genotype)+...etc
///
/// the first term is a product based on error probabilities:
/// P(obs=ACT|AAA) = (1-P(e1))*(P(e2)/3)*(P(e3)/3)
///
/// the second term is a product of genotype base frequencies:
/// P(AAA|genotype=AT) = .5**3
///
/// The likelihood calculated in the function below is an equivelent
/// factorization of the calculation described above.
///
void
pprob_digt_caller::
position_snp_call_pprob_digt(
    const blt_options& opt,
    const extended_pos_info& epi,
    diploid_genotype& dgt,
    const bool is_always_test) const
{
    const snp_pos_info& pi(epi.pi);

    if (pi.get_ref_base()=='N') return;

    dgt.ref_gt=base_to_id(pi.get_ref_base());

    // check that a non-reference call meeting quality criteria even exists:
    if (! is_always_test)
    {
        if (is_spi_allref(pi,dgt.ref_gt)) return;
    }

    // don't spend time on the het bias model for haploid sites:
    const bool useHetVariantFrequencyExtension((! dgt.is_haploid()) && opt.isHetVariantFrequencyExtensionDefined());

    // get likelihood of each genotype
    blt_float_t lhood[DIGT::SIZE];
    get_diploid_gt_lhood(opt,epi,useHetVariantFrequencyExtension,opt.hetVariantFrequencyExtension,lhood);

    // set phredLoghood:
    {
        unsigned gtcount(DIGT::SIZE);
        if (dgt.is_haploid()) gtcount=N_BASE;
        unsigned maxIndex(0);
        for (unsigned gt(1); gt<gtcount; ++gt)
        {
            if (lhood[gt] > lhood[maxIndex]) maxIndex = gt;
        }
        for (unsigned gt(0); gt<gtcount; ++gt)
        {
            dgt.phredLoghood[gt] = ln_error_prob_to_qphred(lhood[gt]-lhood[maxIndex]);
        }
    }


    // get genomic site results:
    calculate_result_set(lhood,lnprior_genomic(dgt.ref_gt,dgt.is_haploid()),dgt.ref_gt,dgt.genome);

    // get polymorphic site results:
    calculate_result_set(lhood,lnprior_polymorphic(dgt.ref_gt,dgt.is_haploid()),dgt.ref_gt,dgt.poly);

    // compute strand-bias here:
    const bool is_compute_sb(true);
    if (is_compute_sb && dgt.is_snp())
    {
        blt_float_t lhood_fwd[DIGT::SIZE];
        get_diploid_gt_lhood(opt,epi,useHetVariantFrequencyExtension,opt.hetVariantFrequencyExtension,lhood_fwd,true,true);
        blt_float_t lhood_rev[DIGT::SIZE];
        get_diploid_gt_lhood(opt,epi,useHetVariantFrequencyExtension,opt.hetVariantFrequencyExtension,lhood_rev,true,false);

        // If max_gt is equal to reference, then go ahead and use it
        // for consistency, even though this makes the SB value
        // useless:
        const unsigned tgt(dgt.genome.max_gt);
        dgt.strand_bias=std::max(lhood_fwd[tgt],lhood_rev[tgt])-lhood[tgt];
    }
    else
    {
        dgt.strand_bias = 0;
    }
}



std::ostream& operator<<(std::ostream& os,const diploid_genotype& dgt)
{
    const result_set& ge(dgt.genome);
    const result_set& po(dgt.poly);

    os << " Q(snp): " << ge.snp_qphred
       << " max_gt: " << DIGT::label(ge.max_gt)
       << " Q(max_gt): "  << ge.max_gt_qphred
       << " max_gt|poly_site: " << DIGT::label(po.max_gt)
       << " Q(max_gt|poly_site): " << po.max_gt_qphred;

    return os;
}
