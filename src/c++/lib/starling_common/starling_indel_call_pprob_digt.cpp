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

/// \file
/// \author Chris Saunders
///

#include "starling_common/AlleleReportInfoUtil.hh"
#include "starling_common/readMappingAdjustmentUtil.hh"
#include "starling_common/starling_indel_call_pprob_digt.hh"

#include <iostream>


//#define DEBUG_INDEL_CALL

#ifdef DEBUG_INDEL_CALL
#include "blt_util/log.hh"
#endif



void
get_het_observed_allele_ratio(
    const unsigned read_length,
    const unsigned min_overlap,
    const IndelKey& indelKey,
    const double het_allele_ratio,
    double& log_ref_prob,
    double& log_indel_prob)
{
    assert(indelKey.type==INDEL::INDEL);

    // the expected relative read depth for two breakpoints separated by a distance of 0:
    const unsigned base_expect( (read_length+1)<(2*min_overlap) ? 0 : (read_length+1)-(2*min_overlap) );

    // Get expected relative read depth for the shorter and longer
    // paths of a general sequence replacement. Note this includes
    // basic insertions and deletions, in these cases
    // spath_break_distance is 0 and spath_expect equals base_expect:
    //
    const double ref_path_expect(base_expect+std::min(indelKey.delete_length(),base_expect));
    const double indel_path_expect(base_expect+std::min(indelKey.insert_length(),base_expect));
    const double ref_path_term((1-het_allele_ratio)*ref_path_expect);
    const double indel_path_term(het_allele_ratio*indel_path_expect);
    const double total_path_term(ref_path_term+indel_path_term);

    if (total_path_term>0)
    {
        const double indel_prob(indel_path_term/total_path_term);
        log_ref_prob=std::log(1.-indel_prob);
        log_indel_prob=std::log(indel_prob);
    }
}



void
get_high_low_het_ratio_lhood(
    const starling_base_options& /*opt*/,
    const starling_base_deriv_options& dopt,
    const starling_sample_options& sample_opt,
    const IndelKey& indelKey,
    const IndelSampleData& indelSampleData,
    const double het_ratio,
    const bool is_tier2_pass,
    const bool is_use_alt_indel,
    double& het_lhood_high,
    double& het_lhood_low)
{
    // handle het ratio and its complement in one step:
    const double chet_ratio(1.-het_ratio);

    const double log_het_ratio(std::log(het_ratio));
    const double log_chet_ratio(std::log(chet_ratio));

    const bool is_breakpoint(indelKey.is_breakpoint());

    het_lhood_high=0;
    het_lhood_low=0;

    for (const auto& score : indelSampleData.read_path_lnp)
    {
        const ReadPathScores& path_lnp(score.second);

        // optionally skip tier2 data:
        if ((! is_tier2_pass) && (! path_lnp.is_tier1_read)) continue;

        // get alt path lnp:
        double alt_path_lnp(path_lnp.ref);
#if 0
        if (is_use_alt_indel && path_lnp.is_alt &&
            (path_lnp.alt > alt_path_lnp))
        {
            alt_path_lnp=path_lnp.alt;
        }
#else
        if (is_use_alt_indel)
        {
            for (const auto& alt : path_lnp.alt_indel)
            {
                if (alt.second>alt_path_lnp) alt_path_lnp=alt.second;
            }
        }
#endif

        const double noindel_lnp(alt_path_lnp);
        const double hom_lnp(path_lnp.indel);

        // allele ratio convention is that the indel occurs at the
        // het_allele ratio and the alternate allele occurs at
        // (1-het_allele_ratio):
        {
            double log_ref_prob(log_chet_ratio);
            double log_indel_prob(log_het_ratio);
            if (! is_breakpoint)
            {
                get_het_observed_allele_ratio(path_lnp.read_length,sample_opt.min_read_bp_flank,
                                              indelKey,het_ratio,log_ref_prob,log_indel_prob);
            }
            const double het_lnp(getLogSum(noindel_lnp+log_ref_prob, hom_lnp+log_indel_prob));

            het_lhood_low += integrateOutMappingStatus(dopt, path_lnp.nonAmbiguousBasesInRead, het_lnp, is_tier2_pass);
        }

        {
            double log_ref_prob(log_het_ratio);
            double log_indel_prob(log_chet_ratio);
            if (! is_breakpoint)
            {
                get_het_observed_allele_ratio(path_lnp.read_length,sample_opt.min_read_bp_flank,
                                              indelKey,chet_ratio,log_ref_prob,log_indel_prob);
            }
            const double het_lnp(getLogSum(noindel_lnp+log_ref_prob, hom_lnp+log_indel_prob));

            het_lhood_high += integrateOutMappingStatus(dopt, path_lnp.nonAmbiguousBasesInRead, het_lnp, is_tier2_pass);
        }
    }
}



static
void
increment_het_ratio_lhood(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const starling_sample_options& sample_opt,
    const IndelKey& indelKey,
    const IndelSampleData& indelSampleData,
    const double het_ratio,
    const bool is_tier2_pass,
    const bool is_use_alt_indel,
    double* const lhood)
{
    // high and low allele ratio variants:
    double het_lhood_high;
    double het_lhood_low;

    get_high_low_het_ratio_lhood(opt,dopt,sample_opt,
                                 indelKey,indelSampleData,het_ratio,is_tier2_pass,
                                 is_use_alt_indel,
                                 het_lhood_high,het_lhood_low);

    lhood[STAR_DIINDEL::HET] = getLogSum(lhood[STAR_DIINDEL::HET], het_lhood_low, het_lhood_high);
}



void
get_sum_path_pprob(
    const starling_base_deriv_options& dopt,
    const IndelSampleData& indelSampleData,
    const bool is_tier2_pass,
    const bool is_use_alt_indel,
    ReadPathScores& total_pprob,
    const bool is_init_total)
{
    static const double initval(0);

    if (is_init_total)
    {
        total_pprob.ref=initval;
        total_pprob.indel=initval;
        total_pprob.nonAmbiguousBasesInRead=0;
    }

    typedef std::map<IndelKey,unsigned> aimap_t;
    aimap_t alt_indel_index;

    for (const auto& score : indelSampleData.read_path_lnp)
    {
        const ReadPathScores& path_lnp(score.second);

        // optionally skip tier2 data:
        if ((! is_tier2_pass) && (! path_lnp.is_tier1_read)) continue;

        const ReadPathScores path_pprob(indel_lnp_to_pprob(dopt,path_lnp,is_tier2_pass,is_use_alt_indel));

        total_pprob.indel += path_pprob.indel;
        total_pprob.ref += path_pprob.ref;

        if (! is_use_alt_indel) continue;

        for (const auto& alt : path_pprob.alt_indel)
        {
            aimap_t::iterator tj(alt_indel_index.find(alt.first));
            if (tj==alt_indel_index.end())
            {
                alt_indel_index[alt.first]=total_pprob.alt_indel.size();
                total_pprob.alt_indel.push_back(alt);
            }
            else
            {
                total_pprob.alt_indel[tj->second].second += alt.second;
            }
        }
    }
}



void
get_indel_digt_lhood(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const starling_sample_options& sample_opt,
    const IndelKey& indelKey,
    const IndelSampleData& indelSampleData,
    const bool is_tier2_pass,
    const bool is_use_alt_indel,
    double* const lhood)
{
    static const double loghalf(-std::log(2.));

    for (unsigned gt(0); gt<STAR_DIINDEL::SIZE; ++gt) lhood[gt] = 0.;

    const bool is_breakpoint(indelKey.is_breakpoint());

    for (const auto& score : indelSampleData.read_path_lnp)
    {
        const ReadPathScores& path_lnp(score.second);

        // optionally skip tier2 data:
        if ((! is_tier2_pass) && (! path_lnp.is_tier1_read)) continue;

        // get alt path lnp:
        double alt_path_lnp(path_lnp.ref);
#if 0
        if (is_use_alt_indel && path_lnp.is_alt &&
            (path_lnp.alt > alt_path_lnp))
        {
            alt_path_lnp=path_lnp.alt;
        }
#else
        if (is_use_alt_indel)
        {
            for (const auto& alt : path_lnp.alt_indel)
            {
                if (alt.second>alt_path_lnp) alt_path_lnp=alt.second;
            }
        }
#endif

        const double noindel_lnp(alt_path_lnp);
        const double hom_lnp(path_lnp.indel);

        // allele ratio convention is that the indel occurs at the
        // het_allele ratio and the alternate allele occurs at
        // (1-het_allele_ratio):

        double log_ref_prob(loghalf);
        double log_indel_prob(loghalf);
        if (! is_breakpoint)
        {
            static const double het_allele_ratio(0.5);
            get_het_observed_allele_ratio(path_lnp.read_length,sample_opt.min_read_bp_flank,
                                          indelKey,het_allele_ratio,log_ref_prob,log_indel_prob);
        }
        const double het_lnp(getLogSum(noindel_lnp+log_ref_prob,hom_lnp+log_indel_prob));

        lhood[STAR_DIINDEL::NOINDEL] += integrateOutMappingStatus(dopt, path_lnp.nonAmbiguousBasesInRead, noindel_lnp,
                                                                  is_tier2_pass);
        lhood[STAR_DIINDEL::HOM]     += integrateOutMappingStatus(dopt, path_lnp.nonAmbiguousBasesInRead, hom_lnp,
                                                                  is_tier2_pass);
        lhood[STAR_DIINDEL::HET]     += integrateOutMappingStatus(dopt, path_lnp.nonAmbiguousBasesInRead, het_lnp,
                                                                  is_tier2_pass);

#ifdef DEBUG_INDEL_CALL
        //log_os << std::setprecision(8);
        //log_os << "INDEL_CALL i,ref_lnp,indel_lnp,lhood(noindel),lhood(hom),lhood(het): " << i << " " << path_lnp.ref << " " << path_lnp.indel << " " << lhood[STAR_DIINDEL::NOINDEL] << " " << lhood[STAR_DIINDEL::HOM] << " " << lhood[STAR_DIINDEL::HET] << "\n";
#endif
    }



    // mothballing het-bias feature as implemented in this indel model, but keeping as a template for
    // implementation in the new diplotype model. het_bias doc is:
    //
    // "                    - Set bias term for the heterozygous state in the bindel model, such that\n"
    // "                      hets are expected at allele ratios in the range [0.5-x,0.5+x] (default: 0)\n"
    static const bool useHetVariantFrequencyExtension(false);
    static const double hetVariantFrequencyExtension(0.);
    if (useHetVariantFrequencyExtension)
    {
        // loop is currently setup to assume a uniform het ratio subgenotype prior
        const unsigned n_bias_steps(1+static_cast<unsigned>(hetVariantFrequencyExtension/opt.maxHetVariantFrequencyIncrement));
        const double ratio_increment(hetVariantFrequencyExtension/static_cast<double>(n_bias_steps));
        for (unsigned step(0); step<n_bias_steps; ++step)
        {
            const double het_ratio(0.5+(step+1)*ratio_increment);
            increment_het_ratio_lhood(opt,dopt,sample_opt,
                                      indelKey,indelSampleData,het_ratio,is_tier2_pass,is_use_alt_indel,lhood);
        }

        const unsigned n_het_subgt(1+2*n_bias_steps);
        const double subgt_log_prior(std::log(static_cast<double>(n_het_subgt)));
        lhood[STAR_DIINDEL::HET] -= subgt_log_prior;
    }
}
