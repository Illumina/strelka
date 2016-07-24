// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
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


#include "gvcfAlleleInfo.hh"
#include "blt_util/math_util.hh"
#include "common/Exceptions.hh"

#include "boost/math/distributions/binomial.hpp"
#include "boost/math/distributions.hpp"
#include "rnaVariantEmpiricalScoringFeatures.hh"

#include <iostream>
#include <map>
#include <sstream>
#include <typeinfo>





void
GermlineIndelSimpleGenotypeInfo::
set_hap_cigar(
    const unsigned lead,
    const unsigned trail)
{
    using namespace ALIGNPATH;

    cigar.clear();
    if (lead)
    {
        cigar.push_back(path_segment(MATCH,lead));
    }
    if (_indelKey.delete_length())
    {
        cigar.push_back(path_segment(DELETE,_indelKey.delete_length()));
    }
    if (_indelKey.insert_length())
    {
        cigar.push_back(path_segment(INSERT,_indelKey.insert_length()));
    }
    if (trail)
    {
        cigar.push_back(path_segment(MATCH,trail));
    }
}



void
GermlineDiploidIndelSimpleGenotypeInfo::
computeEmpiricalScoringFeatures(
    const bool isRNA,
    const bool isUniformDepthExpected,
    const bool isComputeDevelopmentFeatures,
    const double chromDepth,
    const bool isHetalt)
{
    const double filteredLocusDepth(_indelSampleReportInfo.tier1Depth);
    const double locusDepth(_indelSampleReportInfo.mapqTracker.count);
    const double confidentDepth(_indelSampleReportInfo.total_confident_reads());

    const double chromDepthFactor(safeFrac(1,chromDepth));
    const double filteredLocusDepthFactor(safeFrac(1,filteredLocusDepth));
    const double locusDepthFactor(safeFrac(1,locusDepth));
    const double confidentDepthFactor(safeFrac(1,confidentDepth));

    // cdf of binomial prob of seeing no more than the number of 'allele A' reads out of A reads + B reads, given p=0.5
    // cdf of binomial prob of seeing no more than the number of 'allele B' reads out of A reads + B reads, given p=0.5
    double allelebiaslower;
    double allelebiasupper;
    {
        // allele bias metrics
        const double r0(_indelSampleReportInfo.n_confident_ref_reads);
        const double r1(_indelSampleReportInfo.n_confident_indel_reads);
        const double r2(_indelSampleReportInfo.n_confident_alt_reads);

        if (isHetalt)
        {
            allelebiaslower = cdf(boost::math::binomial(r2 + r1, 0.5), r1);
            allelebiasupper = cdf(boost::math::binomial(r2 + r1, 0.5), r2);
        }
        else
        {
            allelebiaslower = cdf(boost::math::binomial(r0 + r1, 0.5), r0);
            allelebiasupper = cdf(boost::math::binomial(r0 + r1, 0.5), r1);
        }
    }

    if (isRNA)
    {
        features.set(RNA_INDEL_SCORING_FEATURES::QUAL, (_dindel.indel_qphred * chromDepthFactor));
        features.set(RNA_INDEL_SCORING_FEATURES::F_GQX, (gqx * chromDepthFactor));
        features.set(RNA_INDEL_SCORING_FEATURES::REFREP1, (_indelReportInfo.ref_repeat_count));
        features.set(RNA_INDEL_SCORING_FEATURES::IDREP1, (_indelReportInfo.indel_repeat_count));
        features.set(RNA_INDEL_SCORING_FEATURES::RULEN1, (_indelReportInfo.repeat_unit.length()));
        features.set(RNA_INDEL_SCORING_FEATURES::AD0,
                     (_indelSampleReportInfo.n_confident_ref_reads * chromDepthFactor));
        features.set(RNA_INDEL_SCORING_FEATURES::AD1,
                     (_indelSampleReportInfo.n_confident_indel_reads * chromDepthFactor));
        features.set(RNA_INDEL_SCORING_FEATURES::AD2,
                     (_indelSampleReportInfo.n_confident_alt_reads * chromDepthFactor));
        features.set(RNA_INDEL_SCORING_FEATURES::F_DPI, (_indelSampleReportInfo.tier1Depth * chromDepthFactor));

        // +1e-30 to avoid log(0) in extreme cases
        features.set(RNA_INDEL_SCORING_FEATURES::ABlower, (-std::log(allelebiaslower + 1.e-30)));
        features.set(RNA_INDEL_SCORING_FEATURES::AB,
                     (-std::log(std::min(1., 2. * std::min(allelebiaslower, allelebiasupper)) + 1.e-30)));


        // compute any experimental features not currently used in production
        if (isComputeDevelopmentFeatures)
        {
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQ, (gqx * chromDepthFactor));

            // how unreliable are the read mappings near this locus?
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_MQ,
                                    (_indelSampleReportInfo.mapqTracker.getRMS()));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::mapqZeroFraction,
                                    (_indelSampleReportInfo.mapqTracker.getZeroFrac()));

            // how noisy is the locus?
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_DPI_NORM,
                                    (filteredLocusDepth * locusDepthFactor));


            // all of the features below are simply renormalized replacements of the current production feature set
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::QUAL_NORM,
                                    (_dindel.indel_qphred * filteredLocusDepthFactor));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQX_NORM,
                                    (gqx * filteredLocusDepthFactor));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQ_NORM,
                                    (gq * filteredLocusDepthFactor));

            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::AD0_NORM,
                                    (_indelSampleReportInfo.n_confident_ref_reads * confidentDepthFactor));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::AD1_NORM,
                                    (_indelSampleReportInfo.n_confident_indel_reads * confidentDepthFactor));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::AD2_NORM,
                                    (_indelSampleReportInfo.n_confident_alt_reads * confidentDepthFactor));

            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::QUAL_EXACT, (_dindel.indel_qphred));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQX_EXACT, (gqx));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQ_EXACT, (gq));
        }
    }
    else
    {
        {
            double genotype(0);
            if (_dindel.max_gt == STAR_DIINDEL::HOM)
            {
                genotype = 1;
            }
            else
            {
                if (_dindel.is_diplotype_model_hetalt) genotype = 2;
            }
            features.set(GERMLINE_INDEL_SCORING_FEATURES::GENO, genotype);
        }

        features.set(GERMLINE_INDEL_SCORING_FEATURES::IDREP1, (_indelReportInfo.indel_repeat_count));
        features.set(GERMLINE_INDEL_SCORING_FEATURES::RULEN1, (_indelReportInfo.repeat_unit.length()));

        // +1e-30 to avoid log(0) in extreme cases
        features.set(GERMLINE_INDEL_SCORING_FEATURES::ABlower, (-std::log(allelebiaslower + 1.e-30)));
        features.set(GERMLINE_INDEL_SCORING_FEATURES::AB,
                     (-std::log(std::min(1., 2. * std::min(allelebiaslower, allelebiasupper)) + 1.e-30)));

        // how unreliable are the read mappings near this locus?
        features.set(GERMLINE_INDEL_SCORING_FEATURES::F_MQ,
                     (_indelSampleReportInfo.mapqTracker.getRMS()));

        features.set(GERMLINE_INDEL_SCORING_FEATURES::AD1_NORM,
                     (_indelSampleReportInfo.n_confident_indel_reads * confidentDepthFactor));

        features.set(GERMLINE_INDEL_SCORING_FEATURES::F_GQX_EXACT, (gqx));

        // compute any experimental features not currently used in production
        if (isComputeDevelopmentFeatures)
        {
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::REFREP1, (_indelReportInfo.ref_repeat_count));

            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::mapqZeroFraction,
                                    (_indelSampleReportInfo.mapqTracker.getZeroFrac()));

            // how noisy is the locus?
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::F_DPI_NORM,
                                    (filteredLocusDepth * locusDepthFactor));

            // how surprising is the depth relative to expect? This is the only value will be modified for exome/targeted runs
            //
            /// TODO: convert this to pvalue based on Poisson distro?
            double relativeLocusDepth(1.);
            if (isUniformDepthExpected)
            {
                relativeLocusDepth = (locusDepth * chromDepthFactor);
            }
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::TDP_NORM, relativeLocusDepth);

            // all of the features below are simply renormalized replacements of the current production feature set
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::QUAL_NORM,
                                    (_dindel.indel_qphred * filteredLocusDepthFactor));
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQX_NORM,
                                    (gqx * filteredLocusDepthFactor));
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQ_NORM,
                                    (gq * filteredLocusDepthFactor));

            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::AD0_NORM,
                                    (_indelSampleReportInfo.n_confident_ref_reads * confidentDepthFactor));
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::AD2_NORM,
                                    (_indelSampleReportInfo.n_confident_alt_reads * confidentDepthFactor));

            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::QUAL_EXACT, (_dindel.indel_qphred));
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQ_EXACT, (gq));
        }
    }
}



std::ostream&
operator<<(std::ostream& os,
           const GermlineVariantSimpleGenotypeInfo& shmod)
{
    os << "gqx: " << shmod.gqx
       << " gq: " << shmod.gq;
    if (typeid(shmod) == typeid(GermlineDiploidIndelSimpleGenotypeInfo))
    {
        auto imod = dynamic_cast<const GermlineDiploidIndelSimpleGenotypeInfo&>(shmod);

        os << " max_gt: " << DIGT::label(imod.max_gt);
    }

    return os;
}

std::ostream&
operator<<(std::ostream& os,
           const GermlineDiploidSiteSimpleGenotypeInfo& smod)
{
    os << static_cast<GermlineVariantSimpleGenotypeInfo>(smod) << '\n';

    os << "is_unknown: " << smod.is_unknown;
    os << " is_covered: " << smod.is_covered;
    os << " is_used_coverage: " << smod.is_used_covered;
    os << " is_zero_ploidy: " << smod.is_zero_ploidy;

    if (smod.modified_gt != MODIFIED_SITE_GT::NONE)
    {
        os << " modgt: " << MODIFIED_SITE_GT::get_label(smod.modified_gt);
    }

    return os;
}



std::ostream&
operator<<(
    std::ostream& os,
    const GermlineIndelSimpleGenotypeInfo& shi)
{
    os << static_cast<GermlineVariantSimpleGenotypeInfo>(shi) << '\n';

    os << "IndelKey: " << shi._indelKey << "\n";
    //os << "indel_data: " << shi._id << "\n";
    os << "indel_report_info: " << shi._indelReportInfo << "\n";
    os << "indel_sample_info: " << shi._indelSampleReportInfo << "\n";
    os << "cigar: " << shi.cigar << "\n";

    return os;
}



std::ostream&
operator<<(
    std::ostream& os,
    const GermlineDiploidIndelSimpleGenotypeInfo& dic)
{
    os << static_cast<GermlineIndelSimpleGenotypeInfo>(dic) << '\n';

    dic._dindel.dump(os);

    os << "EVS: " << dic.empiricalVariantScore << " max_gt: " << dic.max_gt << "\n";

    return os;
}
