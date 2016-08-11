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


#include "gvcf_locus_info.hh"
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
GermlineFilterKeeper::
write(std::ostream& os) const
{
    if (filters.none())
    {
        os << "PASS";
        return;
    }

    bool is_sep(false);
    for (unsigned i(0); i<GERMLINE_VARIANT_VCF_FILTERS::SIZE; ++i)
    {
        if (! filters.test(i)) continue;

        if (is_sep)
        {
            os << ";";
        }
        else
        {
            is_sep=true;
        }
        os << GERMLINE_VARIANT_VCF_FILTERS::get_label(i);
    }
}



pos_t
GermlineDiploidIndelLocusInfo::
end() const
{
    pos_t result = 0;
    for (auto& x : altAlleles)
        result = std::max(result, x._indelKey.right_pos());
    return result;
}



static
void
add_cigar_to_ploidy(
    const ALIGNPATH::path_t& apath,
    std::vector<unsigned>& ploidy)
{
    using namespace ALIGNPATH;
    int offset(-1);
    for (const auto& ps : apath)
    {
        if (is_segment_align_match(ps.type))
        {
            for (unsigned j(0); j<ps.length; ++j)
            {
                if (offset>=0) ploidy[offset]++;
                offset++;
            }
        }
        else if (ps.type==DELETE)
        {
            offset+=ps.length;
        }
    }
}



void
GermlineDiploidIndelLocusInfo::
add_overlap(
    const reference_contig_segment& ref,
    GermlineDiploidIndelLocusInfo& overlap)
{
    assert(altAlleles.size() == 1);
    assert(overlap.altAlleles.size() == 1);

    auto& firstAllele(getFirstAltAllele());
    auto& overlapAllele(overlap.getFirstAltAllele());

    // there's going to be 1 (possibly empty) fill range in front of one haplotype
    // and one possibly empty fill range on the back of one haplotype
    std::string leading_seq,trailing_seq;
    auto indel_end_pos=std::max(overlapAllele._indelKey.right_pos(),firstAllele._indelKey.right_pos());

    const pos_t indel_begin_pos(pos-1);

    sitePloidy.resize(indel_end_pos-pos,0);

    auto munge_indel = [&] (GermlineDiploidIndelLocusInfo& ii)
    {
        auto& this_call(ii.getFirstAltAllele());
        // extend leading sequence start back 1 for vcf compat, and end back 1 to concat with vcf_indel_seq
        ref.get_substring(indel_begin_pos,(ii.pos-indel_begin_pos)-1,leading_seq);
        const unsigned trail_len(indel_end_pos-this_call._indelKey.right_pos());
        ref.get_substring(indel_end_pos-trail_len,trail_len,trailing_seq);

        this_call.set_hap_cigar(leading_seq.size()+1,
                                trailing_seq.size());

        // add to the ploidy object:
        add_cigar_to_ploidy(this_call.cigar,sitePloidy);
    };
    munge_indel(*this);
    munge_indel(overlap);

    // we only combine pairs of simple het indels on different haplotpyes, so this assertion must hold:
    for (const unsigned pl : sitePloidy)
    {
        assert(pl<2);
    }

    //reduce qual and gt to the lowest of the set:
    firstAllele._dindel.indel_qphred = std::min(firstAllele._dindel.indel_qphred, overlap.getFirstAltAllele()._dindel.indel_qphred);
    firstAllele._dindel.max_gt_qphred = std::min(firstAllele._dindel.max_gt_qphred, overlap.getFirstAltAllele()._dindel.max_gt_qphred);


    // combine filter flags from overlapping loci:
    filters.merge(overlap.filters);

    // combine EVS values. Since the "unset" value is -1, this complex logic is necessary
    if (empiricalVariantScore <0)
        empiricalVariantScore = overlap.empiricalVariantScore;
    else if (overlap.empiricalVariantScore >= 0)
        empiricalVariantScore = std::min(empiricalVariantScore, overlap.empiricalVariantScore);
    firstAllele.gqx = std::min(firstAllele.gqx, overlap.getFirstAltAllele().gqx);

    const unsigned sampleCount(getSampleCount());
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        auto& sampleInfo(getSample(sampleIndex));
        sampleInfo.gq = std::min(sampleInfo.gq, overlap.getSample(sampleIndex).gq);
    }
    altAlleles.push_back(overlap.getFirstAltAllele());
}



void
GermlineDiploidIndelLocusInfo::
getPloidyError(
    const unsigned offset) const
{
    using namespace illumina::common;

    std::ostringstream oss;
    oss << "ERROR: get_ploidy offset '" << offset << "' exceeds ploidy region size '" << sitePloidy.size() << "'\n";
    BOOST_THROW_EXCEPTION(LogicException(oss.str()));
}



void
GermlineDiploidSiteLocusInfo::
computeEmpiricalScoringFeatures(
    const bool isRNA,
    const bool isUniformDepthExpected,
    const bool isComputeDevelopmentFeatures,
    const double chromDepth)
{
    const double chromDepthFactor(safeFrac(1, chromDepth));

    ///TODO STREL-125 generalize to multi-sample
    const auto& firstSampleInfo(getSample(0));

    const double filteredLocusDepth(n_used_calls);
    const double locusDepth(mapqCount);

    const double filteredLocusDepthFactor(safeFrac(1, filteredLocusDepth));
    const double locusDepthFactor(safeFrac(1, locusDepth));

    // get the alt base id (choose second in case of an alt het....)
    unsigned altBase(N_BASE);
    for (unsigned b(0); b < N_BASE; ++b)
    {
        if (b == dgt.ref_gt) continue;
        if (DIGT::expect2(b, allele.max_gt))
        {
            altBase = b;
        }
    }
    assert(altBase != N_BASE);

    const unsigned r0 = alleleObservationCounts(dgt.ref_gt);
    const unsigned r1 = alleleObservationCounts(altBase);

    const double mapqZeroFraction(safeFrac(mapqZeroCount, mapqCount));

    const double locusUsedDepthFraction(filteredLocusDepth * locusDepthFactor);


    if (isRNA)
    {
        double genotype(2.0);
        if (is_het() or is_hetalt()) genotype = 1.0;
        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::GT, (genotype));

        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::QUAL, (dgt.genome.snp_qphred * chromDepthFactor));
        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::F_DP, (n_used_calls * chromDepthFactor));
        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::F_DPF, (n_unused_calls * chromDepthFactor));
        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::F_GQ, (firstSampleInfo.gq * chromDepthFactor));
        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::F_GQX, (allele.gqx * chromDepthFactor));

        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::I_AvgBaseQ, (avgBaseQ));
        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::I_AvgPos, (rawPos));

        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::I_BaseQRankSum, (BaseQRankSum));
        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::I_ReadPosRankSum, (ReadPosRankSum));

        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::I_SNVHPOL, (hpol));
        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::I_SNVSB, (allele.strand_bias));

        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::AD0, (r0 * chromDepthFactor));
        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::AD1, (r1 * chromDepthFactor));

        EVSFeatures.set(RNA_SNV_SCORING_FEATURES::ADR, safeFrac(r0, (r0 + r1)));

        // compute any experimental features not currently used in production
        //
        if (isComputeDevelopmentFeatures)
        {
            EVSDevelopmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::I_MQ, (mapqRMS));
            EVSDevelopmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::I_MQRankSum, (MQRankSum));

            EVSDevelopmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::mapqZeroFraction, mapqZeroFraction);

            EVSDevelopmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::F_DP_NORM, locusUsedDepthFraction);

            EVSDevelopmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::QUAL_NORM,
                                          (dgt.genome.snp_qphred * filteredLocusDepthFactor));
            EVSDevelopmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQX_NORM,
                                          (allele.gqx * filteredLocusDepthFactor));
            EVSDevelopmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQ_NORM,
                                          (firstSampleInfo.gq * filteredLocusDepthFactor));

            EVSDevelopmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::AD0_NORM,
                                          (r0 * filteredLocusDepthFactor));
            EVSDevelopmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::AD1_NORM,
                                          (r1 * filteredLocusDepthFactor));

            EVSDevelopmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::QUAL_EXACT,
                                          (dgt.genome.snp_qphred));
            EVSDevelopmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQX_EXACT, (allele.gqx));
            EVSDevelopmentFeatures.set(RNA_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQ_EXACT, (firstSampleInfo.gq));
        }
    }
    else
    {
        {
            double genotype(0);
            if (is_hetalt()) genotype = 2;
            else if (not is_het()) genotype = 1;
            EVSFeatures.set(GERMLINE_SNV_SCORING_FEATURES::GENO, genotype);
        }

        EVSFeatures.set(GERMLINE_SNV_SCORING_FEATURES::I_MQ, (mapqRMS));
        EVSFeatures.set(GERMLINE_SNV_SCORING_FEATURES::I_SNVHPOL, (hpol));
        EVSFeatures.set(GERMLINE_SNV_SCORING_FEATURES::I_SNVSB, (allele.strand_bias));
        EVSFeatures.set(GERMLINE_SNV_SCORING_FEATURES::I_MQRankSum, (MQRankSum));
        EVSFeatures.set(GERMLINE_SNV_SCORING_FEATURES::I_ReadPosRankSum, (ReadPosRankSum));

        // how surprising is the depth relative to expect? This is the only value will be modified for exome/targeted runs
        /// TODO: convert this to pvalue based on Poisson distro?
        double relativeLocusDepth(1.);
        if (isUniformDepthExpected)
        {
            relativeLocusDepth = (locusDepth * chromDepthFactor);
        }

        EVSFeatures.set(GERMLINE_SNV_SCORING_FEATURES::TDP_NORM, relativeLocusDepth);

        // how noisy is the locus?
        EVSFeatures.set(GERMLINE_SNV_SCORING_FEATURES::F_DP_NORM, locusUsedDepthFraction);

        EVSFeatures.set(GERMLINE_SNV_SCORING_FEATURES::F_GQX_EXACT, (allele.gqx));

        // compute any experimental features not currently used in production
        //
        if (isComputeDevelopmentFeatures)
        {

            // BaseQRankSum
            EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::I_BaseQRankSum, (BaseQRankSum));

            // allele bias metrics
            {
                const double allelebiaslower = cdf(boost::math::binomial(r0 + r1, 0.5), r0);
                const double allelebiasupper = cdf(boost::math::binomial(r0 + r1, 0.5), r1);

                // +1e-30 to avoid log(0) in extreme cases
                EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::ABlower, (-log(allelebiaslower + 1.e-30)));
                EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::AB,
                                              (-log(std::min(1., 2. * std::min(allelebiaslower, allelebiasupper)) + 1.e-30)));
            }

            //The average baseQ of the position of alt allele
            EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::I_RawBaseQ, (avgBaseQ));

            //the average position value within a read of alt allele
            EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::I_RawPos, (rawPos));

            // hom unrelable are the read mappings near this locus?
            EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::mapqZeroFraction, mapqZeroFraction);

            // renormalized features intended to replace the corresponding production feature
            //
            EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::QUAL_NORM,
                                          (dgt.genome.snp_qphred * filteredLocusDepthFactor));
            EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQX_NORM,
                                          (allele.gqx * filteredLocusDepthFactor));
            EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQ_NORM,
                                          (firstSampleInfo.gq * filteredLocusDepthFactor));

            EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::AD0_NORM,
                                          (r0 * filteredLocusDepthFactor));

            EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::QUAL_EXACT,
                                          (dgt.genome.snp_qphred));
       
            EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::F_GQ_EXACT, (firstSampleInfo.gq));

            EVSDevelopmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::AD1_NORM,
                                          (r1 * filteredLocusDepthFactor));
        }
    }
}



std::ostream&
operator<<(std::ostream& os,
           const GermlineDiploidSiteLocusInfo& si)
{
    os << "pos: " << (si.pos+1) << " " << si.get_gt();
    return os;
}



void
GermlineDiploidIndelLocusInfo::
computeEmpiricalScoringFeatures(
    const bool isRNA,
    const bool isUniformDepthExpected,
    const bool isComputeDevelopmentFeatures,
    const double chromDepth)
{
    ///TODO STREL-125 generalize to multi-sample
    const auto& firstSampleInfo(getSample(0));
    const auto& firstIndelSampleInfo(getIndelSample(0));

    /// TODO STREL-125 generalize to multiple alts:
    const auto& firstAltAllele(getFirstAltAllele());

    const auto& sampleReportInfo(firstIndelSampleInfo.reportInfo);

    const double filteredLocusDepth(sampleReportInfo.tier1Depth);
    const double locusDepth(sampleReportInfo.mapqTracker.count);
    const double confidentDepth(sampleReportInfo.total_confident_reads());

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
        const double r0(sampleReportInfo.n_confident_ref_reads);
        const double r1(sampleReportInfo.n_confident_indel_reads);
        const double r2(sampleReportInfo.n_confident_alt_reads);

        if (is_hetalt())
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
        features.set(RNA_INDEL_SCORING_FEATURES::QUAL, (firstAltAllele._dindel.indel_qphred * chromDepthFactor));
        features.set(RNA_INDEL_SCORING_FEATURES::F_GQX, (firstAltAllele.gqx * chromDepthFactor));
        features.set(RNA_INDEL_SCORING_FEATURES::REFREP1, (firstAltAllele._indelReportInfo.ref_repeat_count));
        features.set(RNA_INDEL_SCORING_FEATURES::IDREP1, (firstAltAllele._indelReportInfo.indel_repeat_count));
        features.set(RNA_INDEL_SCORING_FEATURES::RULEN1, (firstAltAllele._indelReportInfo.repeat_unit.length()));
        features.set(RNA_INDEL_SCORING_FEATURES::AD0,
                     (sampleReportInfo.n_confident_ref_reads * chromDepthFactor));
        features.set(RNA_INDEL_SCORING_FEATURES::AD1,
                     (sampleReportInfo.n_confident_indel_reads * chromDepthFactor));
        features.set(RNA_INDEL_SCORING_FEATURES::AD2,
                     (sampleReportInfo.n_confident_alt_reads * chromDepthFactor));
        features.set(RNA_INDEL_SCORING_FEATURES::F_DPI, (sampleReportInfo.tier1Depth * chromDepthFactor));

        // +1e-30 to avoid log(0) in extreme cases
        features.set(RNA_INDEL_SCORING_FEATURES::ABlower, (-std::log(allelebiaslower + 1.e-30)));
        features.set(RNA_INDEL_SCORING_FEATURES::AB,
                     (-std::log(std::min(1., 2. * std::min(allelebiaslower, allelebiasupper)) + 1.e-30)));


        // compute any experimental features not currently used in production
        if (isComputeDevelopmentFeatures)
        {
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQ, (firstAltAllele.gqx * chromDepthFactor));

            // how unreliable are the read mappings near this locus?
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_MQ,
                                    (sampleReportInfo.mapqTracker.getRMS()));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::mapqZeroFraction,
                                    (sampleReportInfo.mapqTracker.getZeroFrac()));

            // how noisy is the locus?
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_DPI_NORM,
                                    (filteredLocusDepth * locusDepthFactor));


            // all of the features below are simply renormalized replacements of the current production feature set
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::QUAL_NORM,
                                    (firstAltAllele._dindel.indel_qphred * filteredLocusDepthFactor));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQX_NORM,
                                    (firstAltAllele.gqx * filteredLocusDepthFactor));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQ_NORM,
                                    (firstSampleInfo.gq * filteredLocusDepthFactor));

            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::AD0_NORM,
                                    (sampleReportInfo.n_confident_ref_reads * confidentDepthFactor));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::AD1_NORM,
                                    (sampleReportInfo.n_confident_indel_reads * confidentDepthFactor));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::AD2_NORM,
                                    (sampleReportInfo.n_confident_alt_reads * confidentDepthFactor));

            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::QUAL_EXACT, (firstAltAllele._dindel.indel_qphred));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQX_EXACT, (firstAltAllele.gqx));
            developmentFeatures.set(RNA_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQ_EXACT, (firstSampleInfo.gq));
        }
    }
    else
    {
        {
            double genotype(0);
            if (firstAltAllele._dindel.max_gt == STAR_DIINDEL::HOM)
            {
                genotype = 1;
            }
            else
            {
                if (firstAltAllele._dindel.is_diplotype_model_hetalt) genotype = 2;
            }
            features.set(GERMLINE_INDEL_SCORING_FEATURES::GENO, genotype);
        }

        features.set(GERMLINE_INDEL_SCORING_FEATURES::IDREP1, (firstAltAllele._indelReportInfo.indel_repeat_count));
        features.set(GERMLINE_INDEL_SCORING_FEATURES::RULEN1, (firstAltAllele._indelReportInfo.repeat_unit.length()));

        // +1e-30 to avoid log(0) in extreme cases
        features.set(GERMLINE_INDEL_SCORING_FEATURES::ABlower, (-std::log(allelebiaslower + 1.e-30)));
        features.set(GERMLINE_INDEL_SCORING_FEATURES::AB,
                     (-std::log(std::min(1., 2. * std::min(allelebiaslower, allelebiasupper)) + 1.e-30)));

        // how unreliable are the read mappings near this locus?
        features.set(GERMLINE_INDEL_SCORING_FEATURES::F_MQ,
                     (sampleReportInfo.mapqTracker.getRMS()));

        features.set(GERMLINE_INDEL_SCORING_FEATURES::AD1_NORM,
                     (sampleReportInfo.n_confident_indel_reads * confidentDepthFactor));

        features.set(GERMLINE_INDEL_SCORING_FEATURES::F_GQX_EXACT, (firstAltAllele.gqx));

        // compute any experimental features not currently used in production
        if (isComputeDevelopmentFeatures)
        {
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::REFREP1, (firstAltAllele._indelReportInfo.ref_repeat_count));

            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::mapqZeroFraction,
                                    (sampleReportInfo.mapqTracker.getZeroFrac()));

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
                                    (firstAltAllele._dindel.indel_qphred * filteredLocusDepthFactor));
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQX_NORM,
                                    (firstAltAllele.gqx * filteredLocusDepthFactor));
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQ_NORM,
                                    (firstSampleInfo.gq * filteredLocusDepthFactor));

            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::AD0_NORM,
                                    (sampleReportInfo.n_confident_ref_reads * confidentDepthFactor));
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::AD2_NORM,
                                    (sampleReportInfo.n_confident_alt_reads * confidentDepthFactor));

            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::QUAL_EXACT, (firstAltAllele._dindel.indel_qphred));
            developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQ_EXACT, (firstSampleInfo.gq));
        }
    }
}



void
GermlineDiploidIndelLocusInfo::
dump(std::ostream& os) const
{
    os << "digt_indel_info\n";
    os << "nCalls: " << altAlleles.size() << " isOverlap: " << _is_overlap << "\n";
    os << "ploidy: ";
    for (const unsigned pl : sitePloidy)
    {
        os << " " << pl;
    }
    os << "\n";
    os << "Calls:\n";
    for (const auto& cl : altAlleles)
    {
        os << cl << "\n";
    }
}
